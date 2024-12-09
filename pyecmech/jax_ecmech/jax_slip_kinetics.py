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

            # if abs_rss_crss_frac > self.t_min:
            #     if abs_rss_crss_frac > self.t_max:
            #         shear_dot= shear_dot.at[islip].set(jnp.copysign(jec.GAM_RATIO_OVFFX * self.gamma_0_w, rss[islip]))
            #     else:
            #         temp = jnp.exp(jnp.log(abs_rss_crss_frac)  * self.inv_exp_m1) * self.gamma_0_w
            #         shear_dot = shear_dot.at[islip].set(temp * rss_crss_frac)
        return shear_dot

    def update_hardness(self, hard_state_0, hard_vals, gdot, delta_time, temp_k):
        evol_vals = self.get_evol_vals(gdot)
        
        init_sol = jnp.zeros_like(hard_state_0)
        args = (hard_state_0, evol_vals, delta_time)

        # solver = snls.SNLSTrDlDenseG(self.compute_resid_jacobian, xtolerance=1e-10, ndim=init_sol.shape[0], args=args)
        # solver.delta_control.deltaInit = 1.0
        # status, xs = solver.solve(init_sol)
        # nfev = solver.nfev

        solver = optx.Dogleg(rtol=1e-6, atol=1e-8, norm=optx.two_norm)
        sol = optx.root_find(self.update_hard_resid, solver=solver, y0=init_sol, args=args, throw=False)
        # jax.debug.print("{}", jnp.linalg.norm(self.update_hard_resid(sol.value, args)))
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
        # !! p_dfac is either zero or blows up
        # !IF (pl%p > one) THEN ! no longer allowed
        # !   mts_dfac = zero
        # !ELSE
        # ! blows up, but just set big
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

        # if self.per_slip_system:
        #     ipss = islip
        # else:
        #     ipss = 0

        # if self.per_slip_system:
        #     xn  = self.xn[ipss]
        #     t_min = self.t_min[ipss]
        #     t_max = self.t_max[ipss]
        # else:
        #     xn = self.xn
        #     t_min = self.t_min
        #     t_max = self.t_max

        gin = values[2 + ipss]
        c_t   = values[2 + ipss + self.num_per_slip]
        gamma_w = values[0]
        gamma_r = values[1]

        gathermal, inv_g = jax.lax.cond(
            self.gathermal,
            lambda: (gin, 1.0 / self.tau_a),
            lambda: (self.tau_a, 1.0 / gin)
        )

        # if self.gathermal:
        #     gathermal = gin
        #     inv_g = 1.0 / self.tau_a
        # else:
        #     gathermal = self.tau_a
        #     inv_g = 1.0 / gin

        athermal_0 = jax.lax.cond(
            jnp.abs(tau) < gathermal,
            lambda: 0.0,
            lambda: (jnp.abs(tau) - gathermal) * inv_g
        )

        # if jnp.abs(tau) < gathermal:
        #     athermal_0 = 0.0
        # else:
        #     athermal_0 = (jnp.abs(tau) - gathermal) * inv_g

        # phonon drag related terms first
        drag_exp_arg = (jnp.abs(tau) - gathermal) / self.phonon_drag_stress

        gdot_r = jax.lax.cond(
            drag_exp_arg < jec.DBL_TINY_SQRT,
            lambda: gamma_r * drag_exp_arg,
            lambda: gamma_r * (1.0 - jnp.exp(-drag_exp_arg))
        )

        # if drag_exp_arg < jec.GAM_RATIO_MIN:
        #     return 0.0
        # if drag_exp_arg < jec.DBL_TINY_SQRT:
        #     gdot_r = gamma_r * drag_exp_arg
        # else:
        #     gdot_r = gamma_r * (1.0 - jnp.exp(-drag_exp_arg))


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

        # if mexp_arg > jec.LN_GAM_RATIO_MIN:
        #     # non-vanishing contribution from balancing MTS-like kinetics
        #     gdot_w -= gamma_w * jnp.exp(mexp_arg)

        gdot_w = jax.lax.cond(
            athermal_0 > t_min,
            lambda: gdot_w + (gamma_w * gdot_w_pl_scaling) * jnp.exp(jnp.log(athermal_0) * xn) * athermal_0,
            lambda: gdot_w
        )

        # if athermal_0 > t_min:
        #     # Related to having a smooth transition between the thermal and phonon drag terms
        #     blog = jnp.log(athermal_0) * xn
        #     gdot_w_power_law = (gamma_w * gdot_w_pl_scaling) * jnp.exp(blog) * athermal_0
        #     gdot_w += gdot_w_power_law

        # Combine thermal and phonon drag terms

        # All the ways we could have returned early but jax doesn't allow that :( 
        # # slip rate is effectively 0 due to thermally activated slip kinetics
        # if pexp_arg < jec.LN_GAM_RATIO_MIN:
        #     return 0.0
        # # purely phonon drag limited slip
        # if athermal_0 > t_max:
        #     return gdot_r * jnp.sign(tau)
        # if drag_exp_arg < jec.GAM_RATIO_MIN:
        #     return 0.0
        # gdot = 1.0 / (1.0 / gdot_w + 1.0 / gdot_r) * jnp.sign(tau)
        # return gdot

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

        # solver = snls.SNLSTrDlDenseG(self.compute_resid_jacobian, xtolerance=1e-10, ndim=init_sol.shape[0], args=args)
        # solver.delta_control.deltaInit = 1.0
        # status, xs = solver.solve(init_sol)

        # nfev = solver.nfev
        solver = optx.Dogleg(rtol=1e-6, atol=1e-8, norm=optx.two_norm)
        sol = optx.root_find(self.update_hard_resid, solver=solver, y0=init_sol, args=args, throw=False)
        # jax.debug.print("{}", jnp.linalg.norm(self.update_hard_resid(sol.value, args)))
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

class SlipKineticOrowanD:
    def __init__(
                 self,
                 params
                 ):

        slip_system_geometry_class = params["slip_system_geometry_class"]
        self.gathermal = params["slip_kinetics_gathermal"]

        self.num_slip_systems = slip_system_geometry_class.num_slip_systems
        self.num_hard = 2 * self.num_slip_systems
        self.isotropic = params["slip_kinetics_isotropic"]
        # if per slip system this affects the C1, C2, and berger's magnitude values
        self.per_slip_system = params["slip_kinetics_per_slip_system"]

        #interaction matrix size we'll determine this from the isotropic param
        if self.isotropic:
            self.nIH = 1
        else:
            self.nIH = self.num_slip_systems ** 2

        if self.per_slip_system:
            self.num_per_slip = self.num_slip_systems
        else:
            self.num_per_slip = 1
        # Our ref_slip_rate, CRSS, C1/T, and b*q_m params
        self.num_vals = 1 + self.num_per_slip + 2 * self.num_slip_systems
        # num_per_slip values are C1, C2, and berger's magnitude values ...
        self.num_params = 12 + 4 * self.num_per_slip + self.num_hard + self.nIH
        # num of evol vals are signed scalar mobile dislocation velocity
        self.num_evolve_vals = self.num_slip_systems

        self.shear_mod_ref = params["shear_mod"]
        self.temp_k_ref    = params["temperature_k_ref"]
        # should be per slip system if option set
        self.bergers_magnitude = params["bergers_magnitude"]
        self.lbar_berg = params["lbar_berg"]
        self.slip_gamma_phonon_ref = params["slip_gamma_phonon_ref"]
        self.phonon_drag_stress = params["phonon_drag_stress"]
        self.attempt_frequency = params["attempt_frequency"]
        # should be per slip system if option set
        self.c1 = params["slip_kinetics_c1"]
        self.tau_a = params["slip_kinetics_peirls_barrier"]
        self.p_exponent = params["slip_kinetics_p_exponent"]
        self.q_exponent = params["slip_kinetics_q_exponent"]
        # should be per slip system if option set
        self.c2 = params["slip_kinetics_c2"]
        self.inter_mat = params["slip_kinetics_interaction_matrix"]
        if self.isotropic:
            self.inter_mat = self.inter_mat * jnp.ones((self.num_slip_systems, self.num_slip_systems))

        xm = 1.0 / (2.0 * ((self.c1 / self.temp_k_ref) * self.shear_mod_ref * self.p_exponent * self.q_exponent))

        self.xnn = 1.0 / xm
        self.xn  = np.atleast_1d(self.xnn - 1.0)
        self.t_min = np.atleast_1d(np.power(jec.GAM_RATIO_MIN, xm))
        self.t_max = np.atleast_1d(np.power(jec.GAM_RATIO_OVF, xm))

        # dislocation evolution stuff
        self.c_ann = params["slip_kinetics_c_annihilation"]
        self.d_ann = params["slip_kinetics_d_annihilation"]
        self.c_trap = params["slip_kinetics_c_trap"]
        self.c_mult = params["slip_kinetics_c_multiplication"]
        self.q_mobile = params["slip_kinetics_q_mobile"]
        self.q_total = params["slip_kinetics_q_total"]
        # Should maybe also allow user set this minimum as well
        self.q_minimum = np.min(self.q_mobile) * 1e-4

        self.forest_matrix = np.zeros((self.num_slip_systems, self.num_slip_systems))

        n = np.copy(slip_system_geometry_class.m_vec).T
        s = np.copy(slip_system_geometry_class.s_vec).T

        for i in range(self.num_slip_systems):
            for j in range(self.num_slip_systems):
                self.forest_matrix[i, j] = 0.5 * (np.abs(n[:, i].dot(s[:, j].T)) + np.abs(n[:, i].dot(np.cross(n[:, j], s[:, j]))))

        self.forest_matrix = jnp.asarray(self.forest_matrix)
        # self.forest_matrix = self.inter_mat

    def get_parameters(self, parameters):

        params["slip_kinetics_gathermal"] = self.gathermal
        params["slip_kinetics_isotropic"] = self.isotropic
        params["slip_kinetics_per_slip_system"] = self.per_slip_system
        params["shear_mod"] = self.shear_mod_ref
        params["temperature_k_ref"] = self.temp_k_ref
        params["bergers_magnitude"] = self.bergers_magnitude
        params["lbar_berg"] = self.lbar_berg
        params["slip_gamma_phonon_ref"] = self.slip_gamma_phonon_ref
        params["phonon_drag_stress"] = self.phonon_drag_stress
        params["attempt_frequency"] = self.attempt_frequency
        params["slip_kinetics_c1"] = self.c1
        params["slip_kinetics_peirls_barrier"] = self.tau_a
        params["slip_kinetics_p_exponent"] = self.p_exponent
        params["slip_kinetics_q_exponent"] = self.q_exponent
        params["slip_kinetics_c2"] = self.c2
        params["slip_kinetics_interaction_matrix"] = self.inter_mat
        params["slip_kinetics_c_annihilation"] = self.c_ann
        params["slip_kinetics_d_annihilation"] = self.d_ann
        params["slip_kinetics_c_trap"] = self.c_trap
        params["slip_kinetics_c_multiplication"] = self.c_mult
        params["slip_kinetics_q_mobile"] = self.q_mobile
        params["slip_kinetics_q_total"] = self.q_total
        
        return params

    def get_history_info(self, names, init, plot, state):
        for i in range(self.num_slip_systems):
            name = "hard_state_qM_" + str(i)
            names.append(name)
            init.append(self.q_mobile[i])
            plot.append(True)
            state.append(True)

        for i in range(self.num_slip_systems):
            name = "hard_state_qT_" + str(i)
            names.append(name)
            init.append(self.q_mobile[i])
            plot.append(True)
            state.append(True)

        return (names, init, plot, state) 

    def get_fixed_reference_rate(self, values):
        return values[0]

    def get_values(self, pressure, temp_k, hard_state):
        iend = self.num_slip_systems * 2
        ibeg = self.num_slip_systems
        int_q = np.sqrt(self.inter_mat.dot(hard_state[ibeg:iend]))
        hdnI  = int_q * self.c2

        inv_sqrth = 1.0 / np.sqrt(hard_state[0:ibeg])
        rate = 1.0 / (1.0 / (self.lbar_berg * self.attempt_frequency * inv_sqrth) + 1.0 / (self.slip_gamma_phonon_ref * hard_state[0:ibeg]))

        values = np.zeros(self.num_vals)
        values[0] = np.max(rate)
        ibeg = 1
        iend = 1 + self.num_slip_systems 
        values[ibeg:iend] = hdnI
        ibeg = iend
        iend = ibeg + self.num_slip_systems 
        values[ibeg:iend] = hard_state[0:self.num_slip_systems]
        ibeg = iend
        iend = ibeg + self.num_per_slip
        values[ibeg:iend] = self.c1 / temp_k

        hard_scale = np.mean(hdnI) 
        return (hard_scale, jnp.asarray(values))

    def mts_inner_calc(self, c_e, denom_i, t_frac):
        # !! p_dfac is either zero or blows up
        # !IF (pl%p > one) THEN ! no longer allowed
        # !   mts_dfac = zero
        # !ELSE
        # ! blows up, but just set big
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

        # if self.per_slip_system:
        #     ipss = islip
        # else:
        #     ipss = 0

        # if self.per_slip_system:
        #     xn  = self.xn[ipss]
        #     t_min = self.t_min[ipss]
        #     t_max = self.t_max[ipss]
        # else:
        #     xn = self.xn
        #     t_min = self.t_min
        #     t_max = self.t_max

        gin = values[1 + islip]
        qm  = values[1 + self.num_slip_systems + islip]
        c_t   = values[1 + 2 * self.num_slip_systems + ipss]
        gamma_w = self.lbar_berg * self.attempt_frequency / jnp.sqrt(qm)
        gamma_r = self.slip_gamma_phonon_ref * qm

        gathermal, inv_g = jax.lax.cond(
            self.gathermal,
            lambda: (gin, 1.0 / self.tau_a),
            lambda: (self.tau_a, 1.0 / gin)
        )

        # if self.gathermal:
        #     gathermal = gin
        #     inv_g = 1.0 / self.tau_a
        # else:
        #     gathermal = self.tau_a
        #     inv_g = 1.0 / gin

        athermal_0 = jax.lax.cond(
            jnp.abs(tau) < gathermal,
            lambda: 0.0,
            lambda: (jnp.abs(tau) - gathermal) * inv_g
        )

        # if jnp.abs(tau) < gathermal:
        #     athermal_0 = 0.0
        # else:
        #     athermal_0 = (jnp.abs(tau) - gathermal) * inv_g

        # phonon drag related terms first
        drag_exp_arg = (jnp.abs(tau) - gathermal) / self.phonon_drag_stress

        gdot_r = jax.lax.cond(
            drag_exp_arg < jec.DBL_TINY_SQRT,
            lambda: gamma_r * drag_exp_arg,
            lambda: gamma_r * (1.0 - jnp.exp(-drag_exp_arg))
        )

        # if drag_exp_arg < jec.GAM_RATIO_MIN:
        #     return 0.0
        # if drag_exp_arg < jec.DBL_TINY_SQRT:
        #     gdot_r = gamma_r * drag_exp_arg
        # else:
        #     gdot_r = gamma_r * (1.0 - jnp.exp(-drag_exp_arg))

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

        # if mexp_arg > jec.LN_GAM_RATIO_MIN:
        #     # non-vanishing contribution from balancing MTS-like kinetics
        #     gdot_w -= gamma_w * jnp.exp(mexp_arg)

        gdot_w = jax.lax.cond(
            athermal_0 > t_min,
            lambda: gdot_w + (gamma_w * gdot_w_pl_scaling) * jnp.exp(jnp.log(athermal_0) * xn) * athermal_0,
            lambda: gdot_w
        )

        # if athermal_0 > t_min:
        #     # Related to having a smooth transition between the thermal and phonon drag terms
        #     blog = jnp.log(athermal_0) * xn
        #     gdot_w_power_law = (gamma_w * gdot_w_pl_scaling) * jnp.exp(blog) * athermal_0
        #     gdot_w += gdot_w_power_law

        # Combine thermal and phonon drag terms

        # All the ways we could have returned early but jax doesn't allow that :( 
        # if tau == 0.0:
        #     return 0.0
        # # slip rate is effectively 0 due to thermally activated slip kinetics
        # if pexp_arg < jec.LN_GAM_RATIO_MIN:
        #     return 0.0
        # # purely phonon drag limited slip
        # if athermal_0 > t_max:
        #     return gdot_r * jnp.sign(tau)
        # if drag_exp_arg < jec.GAM_RATIO_MIN:
        #     return 0.0
        # gdot = 1.0 / (1.0 / gdot_w + 1.0 / gdot_r) * jnp.sign(tau)
        # return gdot

        return jnp.select(condlist=[pexp_arg < jec.LN_GAM_RATIO_MIN or drag_exp_arg < jec.GAM_RATIO_MIN or tau == 0.0, athermal_0 > t_max],
                          choicelist=[0.0, gdot_r *jnp.sign(tau)],
                          default = (1.0 / (1.0 / gdot_w + 1.0 / gdot_r) * jnp.sign(tau)))

    def eval_slip_rates(self, rss, values):
        shear_dot = jnp.zeros(self.num_slip_systems)
        for islip in range(self.num_slip_systems):
            shear_dot = shear_dot.at[islip].set(self.calc_slip_rates(rss[islip], values, islip))
        return shear_dot

    def update_hardness(self, hard_state_0, hard_vals, gdot, delta_time, temp_k):

        hard_state_init = jnp.maximum(hard_state_0, self.q_minimum)
        evol_vals = jnp.abs(gdot) /( hard_state_init[0:self.num_slip_systems] * self.bergers_magnitude)

        hard_state_init = jnp.log(hard_state_init)        
        init_sol = jnp.zeros_like(hard_state_init)
        args = (hard_state_init, evol_vals, delta_time)

        # solver = snls.SNLSTrDlDenseG(self.compute_resid_jacobian, xtolerance=1e-10, ndim=init_sol.shape[0], args=args)
        # solver.delta_control.deltaInit = 1.0
        # status, xs = solver.solve(init_sol)

        # nfev = solver.nfev

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
            \dot{q_M} = \dot{q_{mult}} - \dot{q_{trap}} - \dot{q_{ann}}
            \dot{q_T} = \dot{q_{mult}} - \dot{q_{ann}}

            \dot{q_{mult}} = c_{mult} * q_M * \sqrt{q_{for}} * \nu
            \dot{q_{trap}} = c_{trap} * q_M * \sqrt{q_{for}} * \nu
            \dot{q_{ann}}  = c_{ann} * d_{ann} q_M^2 * q_{for} * \nu
            q_{for} = A_{for} q_T
        '''

        nslip = self.num_slip_systems
        ibeg = nslip
        iend = nslip * 2
        
        hard_state_exp = jnp.exp(hard_state)
        # hard_state_exp = hard_state

        forest_dd = self.forest_matrix.dot(hard_state_exp[ibeg:iend])
        sqrt_forest_dis = jnp.sqrt(forest_dd)

        q_mult_dot = self.c_mult * sqrt_forest_dis * hard_state_exp[0:ibeg] * evol_vals
        q_trap_dot = self.c_trap * sqrt_forest_dis * hard_state_exp[0:ibeg] * evol_vals
        q_ann_dot  = self.c_ann * self.d_ann * hard_state_exp[0:ibeg] * hard_state_exp[0:ibeg] * evol_vals

        q_m_dot = q_mult_dot - q_trap_dot - q_ann_dot
        q_t_dot = q_mult_dot - q_ann_dot

        return jnp.hstack([q_m_dot, q_t_dot])  * (1.0 / hard_state_exp)


class SlipKineticBCCMD:
    def __init__(self, params):
        
        #slip_system_geometry_class = params["slip_system_geometry_class"]
        #self.slip_geom_dynamic = slip_system_geometry_class.dynamic
        slip_system_geometry = params["slip_system_geometry"]
        self.slip_geom_dynamic = (slip_system_geometry == "bcc_pencil")
        
        #self.num_slip_systems = slip_system_geometry_class.num_slip_systems
        self.num_slip_systems = params["num_slip_systems"]
        self.num_hard = self.num_slip_systems
        
        #self.nIH = 1
        # Number of parameters the model needs to be instantiated
        self.num_params = 8+4+1
        # Number of slip kinetic related-variables outputted
        self.num_vals = 2 * self.num_slip_systems + 1
        # num of evol vals are signed scalar mobile dislocation velocity
        self.num_evolve_vals = self.num_hard
        
        self.mu = params["shear_mod"]
        self.bmag = params["bergers_magnitude"]
        
        xm = params["slip_kinetics_exp_m"]
        xnn = 1.0 / xm
        self.xn  = xnn - 1.0
        self.t_min = jnp.power(jec.GAM_RATIO_MIN, xm)
        self.t_max = jnp.power(jec.GAM_RATIO_OVF, xm)
        
        self.gam_w0 = params["slip_kinetics_gam_w0"]
        self.tau_p = params["slip_kinetics_tau_p"]
        self.alpha_p = params["slip_kinetics_alpha_p"]
        self.vmax = params["slip_kinetics_vmax"]
        self.tau_drag = params["slip_kinetics_tau_drag"]
        
        # Hardening paraneters
        self.alpha = params["slip_kinetics_alpha"]
        self.k1 = params["slip_kinetics_k1"]
        self.k2 = params["slip_kinetics_k2"]
        self.krelax = params["slip_kinetics_krelax"]
        self.gdot_0 = params["slip_kinetics_gdot_0"]
        self.temp_k0 = params["slip_kinetics_temp_k0"]
        self.ak = params["slip_kinetics_ak"]
        self.hdn_init = params["slip_kinetics_hdn_init"]
        self.hdn_min = 1e-4 * self.hdn_init
    
    # isn't the argument supposed to be params instead of parameters??
    def get_parameters(self, parameters):
        
        params["shear_mod"] = self.mu
        params["bergers_magnitude"] = self.bmag
        params["slip_kinetics_exp_m"] = self.xm
        params["slip_kinetics_gam_w0"] = self.gam_w0
        params["slip_kinetics_tau_p"] = self.tau_p
        params["slip_kinetics_alpha_p"] = self.alpha_p
        params["slip_kinetics_vmax"] = self.vmax
        params["slip_kinetics_tau_drag"] = self.tau_drag
        params["slip_kinetics_alpha"] = self.alpha
        params["slip_kinetics_k1"] = self.k1
        params["slip_kinetics_k2"] = self.k2
        params["slip_kinetics_krelax"] = self.krelax
        params["slip_kinetics_gdot_0"] = self.gdot_0
        params["slip_kinetics_temp_k0"] = self.temp_k0
        params["slip_kinetics_ak"] = self.ak
        params["slip_kinetics_hdn_init"] = self.hdn_init
        
        return params
    
    def get_history_info(self, names, init, plot, state):
        
        for i in range(self.num_slip_systems):
            name = "rho_" + str(i)
            names.append(name)
            init.append(self.hdn_init)
            plot.append(True)
            state.append(True)
        
        return (names, init, plot, state)
        
    def get_fixed_reference_rate(self, values):
        return self.gam_w0
    
    def get_values(self, pressure, temp_k, hard_state):
        
        values = jnp.zeros(self.num_vals)
        crss = jnp.sum(hard_state[0:self.num_slip_systems])
        crss = self.alpha * self.mu * self.bmag * jnp.sqrt(crss)
        mVals = crss
        for i in range(self.num_slip_systems):
            values = values.at[i].set(crss)
            values = values.at[self.num_slip_systems+i].set(hard_state[i])
        values = values.at[self.num_vals].set(temp_k)
        
        return (mVals, jnp.asarray(values))
        
    def calc_slip_rate(self, tau, chi, crss, rho, tK):
        
        def gdot_fun(tau, chi, crss, rho, tK):
            tau_p = self.tau_p / jnp.cos(chi-self.alpha_p)
            t_eff = jnp.maximum(jnp.abs(tau) - tau_p, 0.0)
            g_i = 1.0 / crss
            t_frac = t_eff * g_i
            t_frac = jnp.copysign(t_frac, tau)
            at = jnp.abs(t_frac)
            gam_w = jax.lax.cond(
                self.gam_w0 < 0.0,
                lambda: jnp.abs(self.gam_w0),
                lambda: rho * self.bmag * self.gam_w0
            )
            gmax = rho * self.bmag * self.vmax * (1.0-jnp.exp(-t_eff/self.tau_drag))
            
            def gdot_fun2(at, gam_w, t_frac, gmax):
                blog = self.xn * jnp.log(at)
                gdot = gam_w * jnp.exp(blog) * t_frac
                # Smooth capping to gmax with Lorentz-like factor
                ac = 10.0
                fact = 1.0 / jnp.power(1.0 + jnp.power(jnp.abs(gdot)/gmax, ac), 1.0/ac)
                gdot = fact * gdot
                return gdot
            
            gdot = jnp.select(
                condlist=[
                    (at > self.t_min) & (at > self.t_max),
                    (at > self.t_min) & (at <= self.t_max)
                ],
                choicelist=[
                    jnp.copysign(jec.GAM_RATIO_OVF * gam_w, tau),
                    gdot_fun2(at, gam_w, t_frac, gmax)
                ],
                default=0.0
            )
            
            return gdot
        
        gdot = jax.lax.cond(
            tau == 0.0,
            lambda: 0.0,
            lambda: gdot_fun(tau, chi, crss, rho, tK)
        )
        
        return gdot
    
    def eval_slip_rates(self, rss, values):
        tK = values[-1]
        shear_dot = jnp.zeros(self.num_slip_systems)
        for islip in range(self.num_slip_systems):
            if rss.size == self.num_slip_systems:
                # to debug only, should always provide chi for this model
                taua, chia = rss[islip], 0.0
            else:
                taua, chia = rss[islip], rss[self.num_slip_systems+islip]
            crss, rhoa = values[islip], values[self.num_slip_systems+islip]
            shear_dot = shear_dot.at[islip].set(self.calc_slip_rate(taua, chia, crss, rhoa, tK))
            #print(' ',islip,'taua',taua,'chia',chia*180.0/jnp.pi,'rhoa',rhoa,'gdot',shear_dot[islip])
        return shear_dot
    
    def update_hardness(self, hard_state_0, hard_vals, gdot, delta_time, temp_k):
        hard_state_init = jnp.log(jnp.maximum(hard_state_0, self.hdn_min))
        evol_vals = jnp.hstack([jnp.abs(gdot), jnp.sum(jnp.abs(gdot))])
        init_sol = jnp.zeros_like(hard_state_init)
        args = (hard_state_init, evol_vals, delta_time, hard_vals, temp_k)
        '''
        solver = snls.SNLSTrDlDenseG(self.compute_resid_jacobian, xtolerance=1e-10, ndim=init_sol.shape[0], args=args)
        solver.delta_control.deltaInit = 1.0
        status, xs = solver.solve(init_sol)
        nfev = solver.nfev
        '''
        self.compute_resid_jacobian(init_sol, args=args)
        solver = optx.Dogleg(rtol=1e-6, atol=1e-8, norm=optx.two_norm)
        sol = optx.root_find(self.update_hard_resid, solver=solver, y0=init_sol, args=args, throw=False)
        xs = sol.value
        nfev = sol.stats["num_steps"]

        x_scale = jnp.minimum(hard_state_init, 1.0)
        hard_delta = xs * x_scale
        hard_state = jnp.exp(hard_state_init + hard_delta)
        
        return (nfev, jnp.copy(hard_state))

    def update_hard_resid(self, x, args=()):
        hard_state_0, evol_vals, delta_time, h_vals, temp_k = args
        x_scale = jnp.minimum(hard_state_0, 1.0)
        res_scale = 1.0 / x_scale
        
        hard_state = hard_state_0 + x * x_scale
        hard_state_dot = self.get_hard_state_dot(hard_state, evol_vals, h_vals, temp_k)
        
        residual = (x * x_scale - hard_state_dot * delta_time) * res_scale
        return residual

    def update_hard_jacob(self, x, args=()):
        return jax.jacfwd(self.update_hard_resid, argnums=0)(x, args)

    def compute_resid_jacobian(self, x, args=()):
        residual = self.update_hard_resid(x, args)
        jacob = self.update_hard_jacob(x, args)
        return (residual, jacob)
    
    def get_hard_state_dot(self, hard_state, evol_vals, h_vals, temp_k):
        nslip = self.num_slip_systems

        hexp = jnp.maximum(jnp.exp(hard_state), jec.EPS)
        gamma = evol_vals[-1]

        def k1_func(xi):
            return self.k1 * (1.0 + (self.ak / (jnp.cos(xi - self.alpha_p))))

        def k2_func():
            gamma_ratio = jax.lax.cond(gamma > jec.GAM_RATIO_MIN, lambda: (gamma / self.gdot_0), lambda: jec.GAM_RATIO_MIN)
            return self.k2 * jnp.log(gamma_ratio) * jnp.log(temp_k / self.temp_k0)

        def f_func(abs_gamma_dot):
            A = 100.0
            t = 0.01
            gamma_ratio = jax.lax.cond(gamma > jec.GAM_RATIO_MIN, lambda: (abs_gamma_dot / gamma), lambda: jec.GAM_RATIO_MAX)
            exp_inner = -A * (gamma_ratio - t)
            return 1.0 - ( 1.0 + jnp.exp(exp_inner))

        def k_relax_func(h):
            relax_term = 1.0 - jnp.exp(- (h - self.hdn_min) / self.hdn_min)
            return relax_term

        a_mat = jnp.eye(nslip)
        amat_rho_sqrt = jnp.sqrt(a_mat.dot(hexp))
        sdot = jnp.zeros(nslip)

        for islip in range(nslip):
            k1 = k1_func(h_vals[islip]) * evol_vals[islip]
            k2 = k2_func() * evol_vals[islip]
            fval = f_func(evol_vals[islip]) * self.krelax * k_relax_func(hexp[islip])
            sdot = sdot.at[islip].set(((k1 * amat_rho_sqrt[islip] - k2 * hexp[islip]) - fval * hexp[islip]) / hexp[islip])

        return sdot

if __name__ == "__main__":

    params = {}
    case = "voce_test"
    case = "oro_test"
    case = "mts_test"

    match case:
        case "oro_test":
            dd_init = 1.0e4

            params["num_slip_systems"] = 12
            params["slip_system_geometry_class"] = jslgeo.SlipGeomFCC(params)
            params["slip_kinetics_gathermal"] = False
            params["slip_kinetics_isotropic"] = True
            params["slip_kinetics_per_slip_system"] = False
            params["shear_mod"] = 1.0
            params["temperature_k_ref"] = 300.0
            params["bergers_magnitude"] = 1.0e-4
            params["lbar_berg"] = 10.0 * params["bergers_magnitude"]
            params["slip_gamma_phonon_ref"] = 1.0e3
            params["phonon_drag_stress"] = 0.02
            params["attempt_frequency"] = 1.0e5
            params["slip_kinetics_c1"] = 20000.0
            params["slip_kinetics_peirls_barrier"] = 0.004
            params["slip_kinetics_p_exponent"] = 0.28
            params["slip_kinetics_q_exponent"] = 1.34
            params["slip_kinetics_c2"] = params["shear_mod"] * params["bergers_magnitude"]
            params["slip_kinetics_interaction_matrix"] = 1.0
            params["slip_kinetics_c_annihilation"] = 2.0e-4
            params["slip_kinetics_d_annihilation"] = 6.0 * params["bergers_magnitude"]
            params["slip_kinetics_c_trap"] = 1.0e-3
            params["slip_kinetics_c_multiplication"] = 2.5e-3
            params["slip_kinetics_q_mobile"] = dd_init * np.ones(params["num_slip_systems"])
            params["slip_kinetics_q_total"] = 4.0 * dd_init * np.ones(params["num_slip_systems"])

            params["slip_gamma_phonon_ref"] *= (1.0 / params["slip_kinetics_q_mobile"][0])
            params["attempt_frequency"] *= np.sqrt(params["slip_kinetics_q_mobile"][0])

            skvpl = SlipKineticOrowanD(params)

            delta_time = 0.001

            gdot = np.ones(12)
            gdot[:] *= 0.1 * np.float64(np.r_[0:12])
            gdot[0] = 1.0
            
            hard_state_0 = np.ones(params["num_slip_systems"] * 2)
            hard_state_0[0:params["num_slip_systems"]] = params["slip_kinetics_q_mobile"]
            hard_state_0[params["num_slip_systems"]:(2*params["num_slip_systems"])] = params["slip_kinetics_q_total"]
            hard_vals = np.zeros(params["num_slip_systems"])
            temp_k = 300.0
            nfev, hard_state = skvpl.update_hardness(jnp.asarray(hard_state_0), jnp.asarray(hard_vals), jnp.asarray(gdot), delta_time, temp_k)
            print(hard_state)

            init_tau = 1.0e-2
            pressure = 0.0
            hard_vals, kin_vals = skvpl.get_values(pressure, temp_k, hard_state_0)

            taua = np.ones(params["num_slip_systems"]) * init_tau
            gdots_update = skvpl.eval_slip_rates(taua, kin_vals)

            print(gdots_update)
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
