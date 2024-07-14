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

import jax_ecmech_const as jec
import jax_slip_geom as jslgeo
import jax_snls as snls

import optimistix as optx

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

        self.hard_state_0 = jnp.asarray([self.crss0])

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
        names.append("hard_state_crss")
        init.append(self.hard_state_0[0])
        plot.append(True)
        state.append(True)

        return (names, init, plot, state)

    def get_fixed_reference_rate(self, values):
        return self.gamma_0_w

    def get_values(self, pressure, temp_k, hard_state):
        values = []
        values.append(hard_state[0])
        return (values[0], jnp.asarray(values))

    def eval_slip_rates(self, rss, values):
        def shear_abv_min(rss, abs_rss_crss_frac, rss_crss_frac):
            return jax.lax.cond(
                    abs_rss_crss_frac > self.t_max,
                    lambda: jnp.copysign(jec.GAM_RATIO_OVFFX * self.gamma_0_w, rss),
                    lambda: jnp.exp(jnp.log(abs_rss_crss_frac)  * self.inv_exp_m1) * self.gamma_0_w * rss_crss_frac
                )
        
        shear_dot = jnp.zeros(self.num_slip_systems)

        for islip in range(self.num_slip_systems):
            inv_crss = 1.0 / values[0]
            rss_crss_frac = rss[islip] * inv_crss
            abs_rss_crss_frac = jnp.abs(rss_crss_frac)

            shear_dot = jax.lax.cond(
                abs_rss_crss_frac > self.t_min,
                lambda: shear_dot.at[islip].set(shear_abv_min(rss[islip], abs_rss_crss_frac, rss_crss_frac)),
                lambda: shear_dot
            )

        return shear_dot

    def update_hardness(self, hard_state_0, hard_vals, gdot, delta_time, temp_k):
        evol_vals = self.get_evol_vals(gdot)
        
        init_sol = jnp.zeros_like(hard_state_0)
        args = (hard_state_0, evol_vals, delta_time)

        solver = optx.Dogleg(rtol=1e-6, atol=1e-8, norm=optx.two_norm)
        sol = optx.root_find(self.update_hard_resid, solver=solver, y0=init_sol, args=args, throw=False)

        xs = sol.value
        nfev = sol.stats["num_steps"]

        x_scale = jnp.minimum(hard_state_0, 1.0)
        hard_state = hard_state_0 + xs * x_scale #res.x * x_scale

        return (nfev, jnp.copy(hard_state))

    def get_evol_vals(self, gdot):
        # recompute effective shear rate here versus using a stored value
        abs_shear_rate_sum = jnp.sum(jnp.abs(gdot))

        crss_sat = jax.lax.cond(
            abs_shear_rate_sum > jec.DBL_TINY_SQRT,
            lambda: self.crss_sat * jnp.power((abs_shear_rate_sum / self.gamma_sat_0), self.exp_m_sat),
            lambda: self.crss_sat
        )

        return jnp.asarray([abs_shear_rate_sum, crss_sat])

    def update_hard_resid(self, x, args=()):
        hard_state_0, evol_vals, delta_time = args
        x_scale = jnp.minimum(hard_state_0, 1.0)
        res_scale = 1.0 / x_scale

        hard_state = hard_state_0 + x * x_scale
        hard_state_dot = self.get_hard_state_dot(hard_state, evol_vals)

        residual = (x * x_scale - hard_state_dot * delta_time) * res_scale

        return residual

    def update_hard_jacob(self, x, args=()):
        return jax.jacfwd(self.update_hard_resid, argnums=0)(x, args)

    def compute_resid_jacobian(self, x, args=()):
        residual = self.update_hard_resid(x, args)
        jacob = self.update_hard_jacob(x, args)
        return (residual, jacob)

    def get_hard_state_dot(self, hard_state, evol_vals):
        '''
            \dot{crss} = h_0 * \frac{(crss_sat - crss)}{crss_sat - crss_0}^n' * \Sum^{nslip} | \dot{\gamma}_j |
        '''

        inv_term = jax.lax.cond(
            evol_vals[1] > jnp.atleast_1d(self.crss0)[0],
            lambda: 1.0 / (evol_vals[1] - self.crss0),
            lambda: 0.0
        )

        voce_inner_term = jnp.power((evol_vals[1] - hard_state[0]) * inv_term, self.exp_n1)
        return self.h0 * voce_inner_term * (evol_vals[1] - hard_state[0]) * inv_term * evol_vals[0]

class SlipKineticMTSKocksMecking:
    def __init__(
                 self,
                 params
                 ):

        slip_system_geometry_class = params["slip_system_geometry_class"]
        self.gathermal = params["slip_kinetics_gathermal"]
        self.num_slip_systems = slip_system_geometry_class.num_slip_systems
        self.num_hard = 1
        # if per slip system this affects the C1, C2, and berger's magnitude values
        self.per_slip_system = params["slip_kinetics_per_slip_system"]

        if self.per_slip_system:
            self.num_per_slip = self.num_slip_systems
        else:
            self.num_per_slip = 1
        # Our ref_slip_rate, CRSS, C1/T, and b*q_m params
        self.num_vals = 2 + self.num_per_slip + self.num_per_slip
        # num_per_slip values are C1, C2, and berger's magnitude values ...
        self.num_params = 8 + 3 * self.num_per_slip + 4 + self.num_hard
        # num of evol vals are signed scalar mobile dislocation velocity
        self.num_evolve_vals = 2

        self.shear_mod_ref = params["shear_mod"]
        self.temp_k_ref    = params["temperature_k_ref"]
        # should be per slip system if option set
        self.slip_gamma_phonon_ref = params["slip_gamma_phonon_ref"]
        self.phonon_drag_stress = params["phonon_drag_stress"]
        self.slip_gamma_thermal_ref = params["slip_gamma_thermal_ref"]
        # should be per slip system if option set
        self.c1 = params["slip_kinetics_c1"]
        self.tau_a = params["slip_kinetics_peirls_barrier"]
        self.p_exponent = params["slip_kinetics_p_exponent"]
        self.q_exponent = params["slip_kinetics_q_exponent"]
        # should be per slip system if option set
        self.g0 = params["slip_kinetics_g0_hard"]
        self.s = params["slip_kinetics_s_hard"]

        xm = 1.0 / (2.0 * ((self.c1 / self.temp_k_ref) * self.shear_mod_ref * self.p_exponent * self.q_exponent))

        self.xnn = 1.0 / xm
        self.xn  = np.atleast_1d(self.xnn - 1.0)
        self.t_min = np.atleast_1d(np.power(jec.GAM_RATIO_MIN, xm))
        self.t_max = np.atleast_1d(np.power(jec.GAM_RATIO_OVF, xm))

        # dislocation evolution stuff
        self.k1 = params["slip_kinetics_k1"]
        self.k2_ref = params["slip_kinetics_k2_ref"]
        self.gamma_ref = params["slip_kinetics_gamma_ref"]
        self.n_inv = params["slip_kinetics_n_inv"]
        self.h0 = params["slip_kinetics_dd_ref"] 

        self.h0_min = self.h0 * 1e-4

    def get_parameters(self, parameters):

        params["slip_kinetics_gathermal"] = self.gathermal
        params["slip_kinetics_per_slip_system"] = self.per_slip_system
        params["shear_mod"] = self.shear_mod_ref
        params["temperature_k_ref"] = self.temp_k_ref
        params["slip_gamma_phonon_ref"] = self.slip_gamma_phonon_ref
        params["slip_gamma_thermal_ref"] = self.slip_gamma_thermal_ref
        params["phonon_drag_stress"] = self.phonon_drag_stress
        params["slip_kinetics_c1"] = self.c1
        params["slip_kinetics_peirls_barrier"] = self.tau_a
        params["slip_kinetics_p_exponent"] = self.p_exponent
        params["slip_kinetics_q_exponent"] = self.q_exponent
        params["slip_kinetics_g0_hard"] = self.g0
        params["slip_kinetics_s_hard"] = self.s
        params["slip_kinetics_k1"] = self.k1
        params["slip_kinetics_k2_ref"] = self.k2_ref
        params["slip_kinetics_gamma_ref"] = self.gamma_ref
        params["slip_kinetics_n_inv"] = self.n_inv
        params["slip_kinetics_dd_ref"] = self.h0

        return params

    def get_history_info(self, names, init, plot, state):
        names.append("hard_state_h0")
        init.append(self.h0)
        plot.append(True)
        state.append(True)
        return (names, init, plot, state) 

    def get_fixed_reference_rate(self, values):
        return values[0] + values[1]

    def get_values(self, pressure, temp_k, hard_state):

        values = np.zeros(self.num_vals)

        sqrt_dd = jnp.sqrt(hard_state[0])

        values[0] = self.slip_gamma_thermal_ref / sqrt_dd
        values[1] = self.slip_gamma_phonon_ref * hard_state[0]

        crss = self.g0 + self.s * sqrt_dd

        values[2:(2 + self.num_per_slip)] = crss
        values[(2 + self.num_per_slip): (2 + 2 * self.num_per_slip)] = self.c1 / temp_k

        hd_scale = np.mean(crss)

        return (hd_scale, values)

    def mts_inner_calc(self, c_e, denom_i, t_frac):
        p_func = jax.lax.cond(
            jnp.abs(t_frac) < jec.DBL_TINY_SQRT,
            lambda: 0.0,
            lambda: np.sign(t_frac) * jnp.power(jnp.abs(t_frac), self.p_exponent)
        )

        q_arg = 1.0 - p_func

        pq_fac = jax.lax.cond(
            q_arg < jec.DBL_TINY_SQRT,
            lambda: 0.0,
            lambda: jnp.sign(q_arg) * jnp.power(jnp.abs(q_arg), self.q_exponent)
        )

        return -c_e * pq_fac

    def calc_slip_rates(self, tau, values, islip):
        if tau == 0.0:
            return 0.0
        # slip_rate
        gdot_w_pl_scaling = 10.0

        ipss = jax.lax.cond(
            self.per_slip_system,
            lambda: islip,
            lambda: 0
        )

        xn = self.xn[ipss]
        t_min = self.t_min[ipss]
        t_max = self.t_max[ipss]

        gin = values[2 + ipss]
        c_t   = values[2 + ipss + self.num_per_slip]
        gamma_w = values[0]
        gamma_r = values[1]

        gathermal, inv_g = jax.lax.cond(
            self.gathermal,
            lambda: (gin, 1.0 / self.tau_a),
            lambda: (self.tau_a, 1.0 / gin)
        )

        athermal_0 = jax.lax.cond(
            jnp.abs(tau) < gathermal,
            lambda: 0.0,
            lambda: (jnp.abs(tau) - gathermal) * inv_g
        )

        # phonon drag related terms first
        drag_exp_arg = (jnp.abs(tau) - gathermal) / self.phonon_drag_stress

        gdot_r = jax.lax.cond(
            drag_exp_arg < jec.DBL_TINY_SQRT,
            lambda: gamma_r * drag_exp_arg,
            lambda: gamma_r * (1.0 - jnp.exp(-drag_exp_arg))
        )

        # thermally activated slip kinetic terms next
        c_e = c_t * self.shear_mod_ref
        pt_frac = (jnp.abs(tau) - gathermal) * inv_g
        pexp_arg = self.mts_inner_calc(c_e, inv_g, pt_frac) 

        gdot_w = gamma_w * jnp.exp(pexp_arg) 

        mt_frac = (-jnp.abs(tau) - gathermal) * inv_g
        mexp_arg = self.mts_inner_calc(c_e, inv_g, mt_frac)

        gdot_w = jax.lax.cond(
            mexp_arg > jec.LN_GAM_RATIO_MIN,
            lambda: gdot_w - gamma_w * jnp.exp(mexp_arg),
            lambda: gdot_w
        )

        gdot_w = jax.lax.cond(
            athermal_0 > t_min,
            lambda: gdot_w + (gamma_w * gdot_w_pl_scaling) * jnp.exp(jnp.log(athermal_0) * xn) * athermal_0,
            lambda: gdot_w
        )

        # Combine thermal and phonon drag terms

        return jnp.select(condlist=[pexp_arg < jec.LN_GAM_RATIO_MIN or drag_exp_arg < jec.GAM_RATIO_MIN, athermal_0 > t_max],
                          choicelist=[0.0, gdot_r *jnp.sign(tau)],
                          default = (1.0 / (1.0 / gdot_w + 1.0 / gdot_r) * jnp.sign(tau)))

    def eval_slip_rates(self, rss, values):
        shear_dot = jnp.zeros(self.num_slip_systems)
        for islip in range(self.num_slip_systems):
            shear_dot = shear_dot.at[islip].set(self.calc_slip_rates(rss[islip], values, islip))
        return shear_dot

    def update_hardness(self, hard_state_0, hard_vals, gdot, delta_time, temp_k):

        shear_eff = jnp.sum(jnp.abs(gdot))
        k2 = jax.lax.cond(
            shear_eff > jec.DBL_TINY_SQRT,
            lambda: self.k2_ref * jnp.power((self.gamma_ref / shear_eff), self.n_inv),
            lambda: self.k2_ref
        )

        evol_vals = jnp.asarray([shear_eff, k2])

        hard_state_init = jnp.maximum(hard_state_0, self.h0_min)
        hard_state_init = jnp.log(hard_state_init)
        init_sol = jnp.zeros_like(hard_state_init)
        args = (hard_state_init, evol_vals, delta_time)

        solver = optx.Dogleg(rtol=1e-6, atol=1e-8, norm=optx.two_norm)
        sol = optx.root_find(self.update_hard_resid, solver=solver, y0=init_sol, args=args, throw=False)

        xs = sol.value
        nfev = sol.stats["num_steps"]
        x_scale = jnp.minimum(hard_state_init, 1.0)
        hard_delta = xs * x_scale
        hard_state = jnp.exp(hard_state_init + hard_delta)

        return (nfev, jnp.copy(hard_state))

    def update_hard_resid(self, x, args=()):
        hard_state_0, evol_vals, delta_time = args
        x_scale = jnp.minimum(hard_state_0, 1.0)
        res_scale = 1.0 / x_scale

        hard_state = hard_state_0 + x * x_scale
        hard_state_dot = self.get_hard_state_dot(hard_state, evol_vals)
        residual = (x * x_scale - hard_state_dot * delta_time) * res_scale

        return residual

    def update_hard_jacob(self, x, args=()):
        return jax.jacfwd(self.update_hard_resid, argnums=0)(x, args)

    def compute_resid_jacobian(self, x, args=()):
        residual = self.update_hard_resid(x, args)
        jacob = self.update_hard_jacob(x, args)
        return (residual, jacob)

    def get_hard_state_dot(self, hard_state, evol_vals):
        '''
            In non-log space
            sdot = (k1 * sqrt(h) - k2 * h) * shear_eff 
        '''
        return (self.k1 * jnp.exp(hard_state * -0.5) - evol_vals[1]) * evol_vals[0]

if __name__ == "__main__":

    params = {}
    case = "voce_test"
    case = "mts_test"

    match case:
        case "voce_test":
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
            hard_vals = np.zeros(12)
            temp_k = 300.0

            nfev, hard_state = skvpl.update_hardness(hard_state_0, hard_vals, gdot, delta_time, temp_k)

            print(nfev)
            print(hard_state, hUpdtVal)

        case "mts_test":
            params["num_slip_systems"] = 12
            params["slip_system_geometry_class"] = jslgeo.SlipGeomFCC(params)
            params["slip_kinetics_gathermal"] = False
            params["slip_kinetics_per_slip_system"] = False
            params["shear_mod"] = 1.0
            params["temperature_k_ref"] = 300.0
            params["slip_gamma_phonon_ref"] = 1.0e3
            params["slip_gamma_thermal_ref"] = 20.0
            params["phonon_drag_stress"] = 0.02
            params["slip_kinetics_c1"] = 20000.0
            params["slip_kinetics_peirls_barrier"] = 0.004
            params["slip_kinetics_p_exponent"] = 0.28
            params["slip_kinetics_q_exponent"] = 1.34
            params["slip_kinetics_g0_hard"] = 10.0e-5
            params["slip_kinetics_s_hard"] = 5.0e-5
            params["slip_kinetics_k1"] = 100.0
            params["slip_kinetics_k2_ref"] = 10.0
            params["slip_kinetics_gamma_ref"] = 1.0e-6
            params["slip_kinetics_n_inv"] = 0.05
            params["slip_kinetics_dd_ref"] = 0.25

            hUpdtVal = 0.6633659171982
            delta_time = 1e-1
            gdot = np.zeros(12)
            gdot[0] = 1.0 / 12

            skvpl = SlipKineticMTSKocksMecking(params)
            hard_state_0 = np.ones(1) * params["slip_kinetics_dd_ref"]
            hard_vals = np.zeros(12)
            temp_k = 300.0

            nfev, hard_state = skvpl.update_hardness(hard_state_0, hard_vals, gdot, delta_time, temp_k)

            print(nfev)
            print(hard_state[0], hUpdtVal)

            init_tau = 10.0e-3
            pressure = 0.0
            hard_vals, kin_vals = skvpl.get_values(pressure, temp_k, hard_state_0)

            taua = np.ones(params["num_slip_systems"]) * init_tau
            gdots_update = skvpl.eval_slip_rates(taua, kin_vals)

            gdot_expected = 64.795444829571
            print(gdots_update[0], gdot_expected)

