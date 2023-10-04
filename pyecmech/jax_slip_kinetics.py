#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 08:22:12 2023

@author: carson16
"""

import numpy as np

import jax
import jax.numpy as jnp
import jax.lax.linalg as lax_linalg
from jax import custom_jvp
from functools import partial
from jax import lax
from jax.numpy.linalg import solve
from jax.config import config; config.update("jax_enable_x64", True)

from scipy.optimize import minimize
import scipy.stats as scist
from scipy.optimize import root

import jax_ecmech_util as jeu
import jax_ecmech_const as jec


class SlipKineticVocePowerLaw:
    def __init__(
                 self,
                 params
                 ):
        self.num_hard = 1
        self.nonlinear = params["slip_kin_nonlinear"]
        self.num_params = 3 + 5 + self.num_hard + self.nonlinear
        self.num_vals = 1
        self.num_val_derivs = 1
        self.num_evolve_vals = 2
        self.num_slip_systems = params["num_slip_systems"]

        # Slip kinetics terms
        self.shear_mod = params["shear_mod"]
        self.exp_m = params["slip_kin_exp_m"]
        self.gamma_0_w = params["slip_kin_gamma_0_w"]
        self.inv_exp_m = 1.0 / self.exp_m
        self.inv_exp_m1 = self.inv_exp_m - 1.0

        self.t_min = jnp.power(jec.GAM_RATIO_MIN, self.inv_exp_m)
        self.t_max = jnp.power(jec.GAM_RATIO_OVF, self.inv_exp_m)

        #Voce hardening terms

        self.h0 = params["slip_kin_h0"]
        self.crss0 = params["slip_kin_crss0"]
        self.crss_sat = params["slip_kin_crss_sat"]
        if self.nonlinear:
            self.exp_n = params["slip_kin_voce_exp_n"]
            self.exp_n1 = self.exp_n - 1.0
        else:
            self.exp_n = 1.0
            self.exp_n1 = 0.0

        self.exp_m_sat = params["slip_kin_voce_exp_m_sat"]
        self.gamma_sat_0 = params["slip_kin_voce_gamma_sat_0"]

    def get_parameters(self, parameters):

        params["slip_kin_nonlinear"] = self.nonlinear
        params["shear_mod"] = self.shear_mod 
        params["slip_kin_exp_m"] = self.exp_m
        params["slip_kin_gamma_0_w"] = self.gamma_0_w
        params["slip_kin_h0"] = self.h0 
        params["slip_kin_crss0"] = self.crss0
        params["slip_kin_crss_sat"] = self.crss_sat
        params["slip_kin_voce_exp_n"] = self.exp_n
        params["slip_kin_voce_exp_m_sat"] = self.exp_m_sat
        params["slip_kin_voce_gamma_sat_0"] = self.gamma_sat_0

        return params

    def get_history_info(self, names, init, plot, state):
        names.append("h")
        init.append(self.crss0)
        plot.append(True)
        state.append(True)

        return (names, init, plot, state) 

    def get_fixed_reference_rate(self, values):
        return self.gamma_0_w

    def get_values(self, pressure, temp_k, hard_state):
        values = []
        values = self.hard_state[0]
        return (values[0], jnp.asarray(values))

    def eval_slip_rates(self, rss, values):
        shear_dot = jnp.zeros(self.num_slip_systems)

        for islip in range(self.num_slip_systems):
            inv_crss = 1.0 / values[0]
            rss_crss_frac = rss[islip] * inv_crss
            abs_rss_crss_frac = jnp.abs(rss_crss_frac)
            if abs_rss_crss_frac > self.t_min:
                if abs_rss_crss_frac > self.t_max:
                    shear_dot.at[islip].set(jnp.copysign(jec.GAM_RATIO_OVFFX * self.gamma_0_w, self.rss[islip]))
                else:
                    temp = jnp.exp(jnp.log(abs_rss_crss_frac)  * self.inv_exp_m1) * self.gamma_0_w
                    shear_dot.at[islip].set(temp * rss_crss_frac)
        return shear_dot


    def update_hardness(self, hard_state_0, hard_vals, gdot, delta_time, temp_k):
        evol_vals = self.get_evol_vals(gdot)
        
        init_sol = jnp.zeros_like(hard_state_0)
        args = (hard_state_0, evol_vals, delta_time)
        res = root(self.update_hard_resid, init_sol, args=args, jac=self.update_hard_jacob, method='hybr', tol=1e-8)

        x_scale = jnp.minimum(hard_state_0, 1.0)
        hard_state = hard_state_0 + res.x * x_scale

        return (res.nfev, jnp.copy(hard_state))

    def get_evol_vals(self, gdot):
        # recompute effective shear rate here versus using a stored value
        abs_shear_rate_sum = jnp.sum(jnp.abs(gdot))

        crss_sat = self.crss_sat
        if abs_shear_rate_sum > jec.DBL_TINY_SQRT:
            crss_sat *= jnp.power((abs_shear_rate_sum / self.gamma_sat_0), self.exp_m_sat)

        return jnp.asarray([abs_shear_rate_sum, crss_sat])

    def update_hard_resid(self, x, hard_state_0, evol_vals, delta_time):
        x_scale = jnp.minimum(hard_state_0, 1.0)
        res_scale = 1.0 / x_scale

        hard_state = hard_state_0 + x * x_scale
        hard_state_dot = self.get_hard_state_dot(hard_state, evol_vals)

        residual = (x * x_scale - hard_state_dot * delta_time) * res_scale

        return residual

    def update_hard_jacob(self, x, hard_state_0, evol_vals, delta_time):
        return jax.jacfwd(self.update_hard_resid, argnums=1)(x, hard_state_0, evol_vals, delta_time)

    def get_hard_state_dot(self, hard_state, evol_vals):
        '''
            \dot{crss} = h_0 * \frac{(crss_sat - crss)}{crss_sat - crss_0}^n' * \Sum^{nslip} | \dot{\gamma}_j |
        '''
        inv_term = 0.0
        if evol_vals[1] > jnp.atleast_1d(self.crss0):
            inv_term = 1.0 / (evol_vals[1] - self.crss0)
        voce_inner_term = jnp.power((evol_vals[1] - hard_state[0]) * inv_term, self.exp_n1)
        return self.h0 * voce_inner_term * (evol_vals[1] - hard_state[0]) * inv_term * evol_vals[0]

if __name__ == "__main__":

    params = {}

    params["slip_kin_nonlinear"] = False
    params["num_slip_systems"] = 12
    params["shear_mod"] = 1.0 
    params["slip_kin_exp_m"] = 0.01
    params["slip_kin_gamma_0_w"] = 1.0
    params["slip_kin_h0"] = 200e-5 
    params["slip_kin_crss0"] = 100e-5
    params["slip_kin_crss_sat"] = 400e-5
    params["slip_kin_voce_exp_n"] = 1.0
    params["slip_kin_voce_exp_m_sat"] = 0.05
    params["slip_kin_voce_gamma_sat_0"] = 1.0e-6

    hUpdtVal_nl = 0.001016575445448
    hUpdtVal = 0.001016620868315
    delta_time = 1e-1
    gdot = np.zeros(12)
    gdot[0] = 1.0 / 12

    skvpl = SlipKineticVocePowerLaw(params)
    hard_state_0 = np.ones(1) * params["slip_kin_crss0"]
    hard_vals = np.zeros(1)
    temp_k = 300.0

    nfev, hard_state = skvpl.update_hardness(hard_state_0, hard_vals, gdot, delta_time, temp_k)

    print(nfev)

    print(hard_state, hUpdtVal)

