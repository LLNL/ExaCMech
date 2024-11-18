#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 08:22:12 2023

@author: carson16
"""

import numpy as np

import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

import jax_ecmech_util as jeu
import jax_ecmech_const as jec

class SlipGeomBase:
    def __init__(self,
                 params):
        self.num_slip_systems = params["num_slip_systems"]
        self.m_vec = jnp.zeros((3, self.num_slip_systems))
        self.s_vec = jnp.zeros((3, self.num_slip_systems))
        self.p_vec = jnp.zeros((jec.NTVEC, self.num_slip_systems))
        self.q_vec = jnp.zeros((jec.NWVEC, self.num_slip_systems))
    
    def get_PQ_chia(self, kirchoff_dev):
        chia = jnp.zeros(self.num_slip_systems)
        return (chia, self.p_vec, self.q_vec)

    def evaluate_RSS(self, kirchoff_dev):
        return jnp.dot(kirchoff_dev[0:-1], self.p_vec)

    def fill_from_mvec_svec(self, mvecs, svecs):
        p_vec = jnp.zeros((jec.NTVEC, self.num_slip_systems))
        q_vec = jnp.zeros((jec.NWVEC, self.num_slip_systems))
        for islip in range(self.num_slip_systems):
            schmid = jnp.outer(svecs[islip, :], mvecs[islip, :])
            pt, qt = jeu.matrix_to_p_q(schmid)
            p_vec = p_vec.at[:, islip].set(pt)
            q_vec = q_vec.at[:, islip].set(qt)
        return (p_vec, q_vec)

    def get_parameters(self, params):
        return params

class SlipGeomFCC(SlipGeomBase):
    def __init__(self,
                 params):
        SlipGeomBase.__init__(self, params)

        self.dynamic = False
        self.num_params = 0
        self.num_slip_systems = 12

        P3 = 1.0 / jnp.sqrt(3.0)
        M3 = -P3
        P2 = 1.0 / jnp.sqrt(2.0)
        M2 = -P2
        Z  = 0.0

        #Slip plane normal CUB111
        mvecs = jnp.asarray([
            [P3, P3, P3],
            [P3, P3, P3],
            [P3, P3, P3],
            [P3, P3, M3],
            [P3, P3, M3],
            [P3, P3, M3],
            [P3, M3, P3],
            [P3, M3, P3],
            [P3, M3, P3],
            [P3, M3, M3],
            [P3, M3, M3],
            [P3, M3, M3]
        ]) 

        svecs = jnp.asarray([
            [Z,  P2, M2],
            [P2, Z,  M2],
            [P2, M2, Z],
            [Z,  P2, P2],
            [P2, Z,  P2],
            [P2, M2, Z],
            [Z,  P2, P2],
            [P2, Z,  M2],
            [P2, P2, Z],
            [Z,  P2, M2],
            [P2, Z,  P2],
            [P2, P2, Z]    
        ])

        self.m_vec = np.copy(mvecs)
        self.s_vec = np.copy(svecs)
        self.p_vec, self.q_vec = self.fill_from_mvec_svec(mvecs, svecs)

class SlipGeomBCC(SlipGeomBase):
    def __init__(self, params):
        SlipGeomBase.__init__(self, params)

        self.dynamic = False
        self.num_params = 1
        self.bcc_type = params["bcc_type"]
        if self.bcc_type == "bcc12":
            self.num_slip_systems = 12
        elif self.bcc_type == "bcc24":
            self.num_slip_systems = 24
        elif self.bcc_type == "bcc48":
            self.num_slip_systems = 48
        else:
            print("Provided invalid bcc_type reverting to bcc12 case")
            self.num_slip_systems = 12
        

        P3 = 1.0 / jnp.sqrt(3.0)
        M3 = -P3
        P2 = 1.0 / jnp.sqrt(2.0)
        M2 = -P2
        Z  = 0.0

        #Slip direction 111
        svecs = jnp.asarray([
            [P3, P3, P3],
            [P3, P3, P3],
            [P3, P3, P3],
            [P3, P3, M3],
            [P3, P3, M3],
            [P3, P3, M3],
            [P3, M3, P3],
            [P3, M3, P3],
            [P3, M3, P3],
            [P3, M3, M3],
            [P3, M3, M3],
            [P3, M3, M3]
        ]) 

        # slip plane normal 110
        mvecs = jnp.asarray([
            [Z,  P2, M2],
            [P2, Z,  M2],
            [P2, M2, Z],
            [Z,  P2, P2],
            [P2, Z,  P2],
            [P2, M2, Z],
            [Z,  P2, P2],
            [P2, Z,  M2],
            [P2, P2, Z],
            [Z,  P2, M2],
            [P2, Z,  P2],
            [P2, P2, Z]    
        ])

        if self.num_slip_systems >= 24:
            P62 = 2.0 / jnp.sqrt(6.0)
            P6  = 1.0 / jnp.sqrt(6.0)
            M62 = -P62
            M6  = -P6

            # slip plane normal 112
            mvecthis = jnp.asarray([
                [M62, P6, P6],
                [P6, M62, P6],
                [P6, P6, M62],
                [M6, M62, P6],
                [P62, P6, P6],
                [M6, P6, M62],
                [P62, M6, P6],
                [M6, P62, P6],
                [M6, M6, M62],
                [P6, P62, P6],
                [M62, M6, P6],
                [P6, M6, M62],                
            ])

            # slip direction 111
            svecthis = jnp.asarray([
                [P3, P3, P3],
                [P3, P3, P3],
                [P3, P3, P3],
                [M3, P3, P3],
                [M3, P3, P3],
                [M3, P3, P3],
                [M3, M3, P3],
                [M3, M3, P3],
                [M3, M3, P3],
                [P3, M3, P3],
                [P3, M3, P3],
                [P3, M3, P3],
            ])

            mvecs = jnp.concatenate((mvecs, mvecthis), axis=0)
            svecs = jnp.concatenate((svecs, svecthis), axis=0)

        if self.num_slip_systems >= 48:

            mPg2a = 1.0 / jnp.sqrt(14.0)
            mPg2b = 2.0 / jnp.sqrt(14.0)
            mPg2c = 3.0 / jnp.sqrt(14.0)

            # 24 {123}<111> slip systems
            mvecthis = jnp.asarray([
                [mPg2c, -mPg2a, -mPg2b],
                [-mPg2b, mPg2c, -mPg2a],
                [-mPg2a, -mPg2b, mPg2c],
                [mPg2a, mPg2c, -mPg2b],
                [-mPg2c, -mPg2b, -mPg2a],
                [mPg2b, -mPg2a, mPg2c],
                [-mPg2c, mPg2a, -mPg2b],
                [mPg2b, -mPg2c, -mPg2a],
                [mPg2a, mPg2b, mPg2c],
                [-mPg2a, -mPg2c, -mPg2b],
                [mPg2c, mPg2b, -mPg2a],
                [-mPg2b, mPg2a, mPg2c],
                [-mPg2a, mPg2c, mPg2b],
                [mPg2c, -mPg2b, mPg2a],
                [-mPg2b, -mPg2a, -mPg2c],
                [-mPg2c, -mPg2a, mPg2b],
                [mPg2b, mPg2c, mPg2a],
                [mPg2a, -mPg2b, -mPg2c],
                [mPg2a, -mPg2c, mPg2b],
                [-mPg2c, mPg2b, mPg2a],
                [mPg2b, mPg2a, -mPg2c],
                [mPg2c, mPg2a, mPg2b],
                [-mPg2b, -mPg2c, mPg2a],
                [-mPg2a, mPg2b, -mPg2c]
            ])

            svecthis = jnp.asarray([
                [P3, P3, P3],
                [P3, P3, P3],
                [P3, P3, P3],
                [M3, P3, P3],
                [M3, P3, P3],
                [M3, P3, P3],
                [M3, M3, P3],
                [M3, M3, P3],
                [M3, M3, P3],
                [P3, M3, P3],
                [P3, M3, P3],
                [P3, M3, P3],
                [P3, P3, M3],
                [P3, P3, M3],
                [P3, P3, M3],
                [M3, P3, M3],
                [M3, P3, M3],
                [M3, P3, M3],
                [M3, M3, M3],
                [M3, M3, M3],
                [M3, M3, M3],
                [P3, M3, M3],
                [P3, M3, M3],
                [P3, M3, M3]
            ])

            mvecs = jnp.concatenate((mvecs, mvecthis), axis=0)
            svecs = jnp.concatenate((svecs, svecthis), axis=0)

        self.m_vec = np.copy(mvecs)
        self.s_vec = np.copy(svecs)

        self.p_vec, self.q_vec = self.fill_from_mvec_svec(mvecs, svecs)
    
    def get_parameters(self, params):
        params["bcc_type"] = self.bcc_type
        return params

if __name__ == "__main__":

    params = {}
    params["num_slip_systems"] = 12
    params["bcc_type"] = "bcc48"

    sgbcc = SlipGeomBCC(params)
