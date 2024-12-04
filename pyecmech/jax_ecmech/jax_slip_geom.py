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
    
    def get_PQ_chia(self, kirchoff_dev, setvals=False):
        chia = jnp.zeros(self.num_slip_systems)
        return (chia, self.p_vec, self.q_vec)

    def evaluate_RSS(self, kirchoff_dev, p_vec):
        return jnp.dot(kirchoff_dev[0:-1], p_vec)

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

        self.dynamic = 0
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

        self.dynamic = 0
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

class SlipGeomBCCPencil(SlipGeomBase):
    def __init__(self, params):
        SlipGeomBase.__init__(self, params)

        self.dynamic = 2 # 0=fixed, 1=fully dynamic, 2=precompute
        self.num_params = 0
        self.num_slip_systems = 4
        
        P3 = 1.0 / jnp.sqrt(3.0)
        M3 = -P3
        P2 = 1.0 / jnp.sqrt(2.0)
        M2 = -P2
        Z  = 0.0
        
        svecs = jnp.asarray([
            [P3, P3, P3],
            [M3, P3, P3],
            [P3, M3, P3],
            [P3, P3, M3],
        ])
        
        mvecs = jnp.asarray([
            [Z, M2, P2],
            [Z, M2, P2],
            [Z, P2, P2],
            [Z, P2, P2],
        ])
        
        self.m_vec = jnp.copy(mvecs)
        self.s_vec = jnp.copy(svecs)
        self.p_vec, self.q_vec = self.fill_from_mvec_svec(mvecs, svecs)
        self.chia = jnp.zeros(self.num_slip_systems)
    
    def _get_PQ_chia(self, kirchoff_dev):
        
        eps = 1e-10
        chia = jnp.zeros(self.num_slip_systems)
        m_vecs = jnp.zeros((self.num_slip_systems,3))
        
        SvecP = jeu.dev_vec_to_sym_vec(kirchoff_dev)
        p = 1.0/3.0 * SvecP[6]
        S = jnp.asarray([[SvecP[0]+p, SvecP[5], SvecP[4]],
                         [SvecP[5], SvecP[1]+p, SvecP[3]],
                         [SvecP[4], SvecP[3], SvecP[2]+p]])
        
        for islip in range(self.num_slip_systems):
            # PK force direction
            #fpk = jnp.zeros(3)
            Sb = jnp.dot(S, self.s_vec[islip,:])
            #if jnp.linalg.norm(Sb) > eps:
            #    fpk = jnp.cross(Sb, self.s_vec[islip,:])
                
            fpk = jax.lax.cond(
                jnp.linalg.norm(Sb) > eps,
                lambda: jnp.cross(Sb, self.s_vec[islip,:]),
                lambda: jnp.zeros(3)
            )
                
            # Normal direction
            '''
            if jnp.linalg.norm(Sb) > eps:
                mvec = jnp.cross(self.s_vec[islip,:], fpk)
                mvec = mvec / jnp.linalg.norm(mvec)
            else:
                mvec = self.m_vec[islip,:]
            '''
            mvec = jax.lax.cond(
                jnp.linalg.norm(Sb) > eps,
                lambda: jnp.cross(self.s_vec[islip,:], fpk),
                lambda: self.m_vec[islip,:]
            )
            mvec = mvec / jnp.linalg.norm(mvec)
            
            m_vecs = m_vecs.at[islip,:].set(mvec)
                
            # MRSSP angle
            n0_vec = 1.0 / jnp.sqrt(2.0) * jnp.array([2.0, -1.0, -1.0]) * self.s_vec[islip,:]
            t0_vec = jnp.cross(n0_vec, self.s_vec[islip,:])
            fx, fy = jnp.dot(fpk, t0_vec), jnp.dot(fpk, n0_vec)
            chi = jnp.arctan2(fy, fx)-jnp.pi/6.0
            
            # Fold into T/AT primary region (-30:30)
            '''
            if chi > 1.0*jnp.pi/6.0 and chi <= 3.0*jnp.pi/6.0:
                chi = jnp.pi/3.0-chi
            elif chi > 3.0*jnp.pi/6.0 and chi <= 5.0*jnp.pi/6.0:
                chi -= 2.0*jnp.pi/3.0
            elif chi >= -7.0*jnp.pi/6.0 and chi < -5.0*jnp.pi/6.0:
                chi = -jnp.pi-chi
            elif chi >= -5.0*jnp.pi/6.0 and chi < -3.0*jnp.pi/6.0:
                chi += 2.0*jnp.pi/3.0
            elif chi >= -3.0*jnp.pi/6.0 and chi < -1.0*jnp.pi/6.0:
                chi = -jnp.pi/3.0-chi
            '''
            chi = jnp.select(
                condlist=[
                    (chi > 1.0*jnp.pi/6.0) & (chi <= 3.0*jnp.pi/6.0),
                    (chi > 3.0*jnp.pi/6.0) & (chi <= 5.0*jnp.pi/6.0),
                    (chi >= -7.0*jnp.pi/6.0) & (chi < -5.0*jnp.pi/6.0),
                    (chi >= -5.0*jnp.pi/6.0) & (chi < -3.0*jnp.pi/6.0),
                    (chi >= -3.0*jnp.pi/6.0) & (chi < -1.0*jnp.pi/6.0)
                ],
                choicelist=[
                    jnp.pi/3.0-chi,
                    chi-2.0*jnp.pi/3.0,
                    -jnp.pi-chi,
                    chi+2.0*jnp.pi/3.0,
                    -jnp.pi/3.0-chi
                ],
                default=chi
            )
            
            chia = chia.at[islip].set(chi)
            
        p_vec, q_vec = self.fill_from_mvec_svec(m_vecs, self.s_vec)
        
        return (chia, p_vec, q_vec)
    
    def get_PQ_chia(self, kirchoff_dev, setvals=False):
        
        if self.dynamic == 2 and not setvals:
            # Return stored values
            chia, p_vec, q_vec = self.chia, self.p_vec, self.q_vec
        else:
            chia, p_vec, q_vec = self._get_PQ_chia(kirchoff_dev)
            
        if setvals:
            # Store values
            self.chia, self.p_vec, self.q_vec = chia, p_vec, q_vec
        
        return (chia, p_vec, q_vec)
    
    

if __name__ == "__main__":

    params = {}
    params["num_slip_systems"] = 12
    params["bcc_type"] = "bcc48"

    sgbcc = SlipGeomBCC(params)
