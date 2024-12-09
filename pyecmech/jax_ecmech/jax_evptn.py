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
jax.config.update("jax_enable_x64", True)

# from scipy.optimize import minimize
# import scipy.stats as scist
# from scipy.optimize import root

import optimistix as optx

import jax_ecmech_util as jeu
import jax_ecmech_const as jec

import jax_eos as jeos
import jax_slip_geom as jslgeo
import jax_slip_kinetics as jslkin
import jax_thermo_elastn as jtelas
import jax_snls as snls

class evptnClass:
    def __init__(
                 self,
                 slip_geom_class,
                 slip_kinetics_class,
                 thermo_elas_class,
                 delta_time,
                 det_vol, 
                 eVref, 
                 pressure_eos,
                 temp_k,
                 hard_state,
                 elas_dev_vec_n,
                 crystal_quat_n,
                 def_dev_vec_samp,
                 spin_dev_vec_samp
                 ):
    
        self.slip_geom_class = slip_geom_class
        self.slip_kinetics_class = slip_kinetics_class
        self.thermo_elas_class = thermo_elas_class
        self.delta_time = delta_time
        self.det_vol = det_vol 
        self.eVref = eVref
        self.pressure_eos = pressure_eos
        self.temp_k = temp_k
        self.hard_state = hard_state 
        self.elas_dev_vec_n = elas_dev_vec_n
        self.crystal_quat_n = crystal_quat_n
        self.def_dev_vec_samp = def_dev_vec_samp 
        self.spin_dev_vec_samp = spin_dev_vec_samp

        self.ind_e_beg = 0
        self.ind_e_end = 5
        self.ind_r_beg = self.ind_e_end
        self.ind_r_end = self.ind_r_beg + 3
        
        self.inv_delta_time = 1.0 / delta_time
        self.inv_det_vol = 1.0 / det_vol
        self.a_vol = jnp.power(det_vol, 1.0 / 3.0)
        self.inv_a_vol = 1.0 / self.a_vol
        
        self.hard_scale, self.kinetic_vals = self.slip_kinetics_class.get_values(self.pressure_eos, self.temp_k, self.hard_state)

        self.elas_scale = jec.ELAS_SCALE
        self.rot_scale  = jec.ROT_SCALE

        adot_ref = self.slip_kinetics_class.get_fixed_reference_rate(self.kinetic_vals)
        eff = jnp.linalg.norm(def_dev_vec_samp)
        
        self.epsdot_scale_inv = jax.lax.cond(
            eff < (jnp.sqrt(jnp.finfo(jnp.float64).eps) * adot_ref),
            lambda: 1.0 / adot_ref,
            lambda: jnp.minimum(1.0 / eff, 1.0e6 * delta_time)
        )
        
        self.rotation_incr_scale_inv = self.inv_delta_time * self.epsdot_scale_inv


    def get_state_from_x(self, x):
        delta_elas_vec = self.elas_scale * x[self.ind_e_beg:self.ind_e_end]
        elas_vec = self.elas_dev_vec_n + delta_elas_vec
        
        exp_map = self.rot_scale * x[self.ind_r_beg:self.ind_r_end]
        dquat = jeu.exp_map_to_quat(exp_map)
        crystal_quat = jeu.update_quat_rot(dquat, self.crystal_quat_n)
        
        return (elas_vec, crystal_quat)
    
    def elas_strain_to_kirchoff_stress(self, elas_dev_vec):
        
        elas_dev = self.inv_a_vol * elas_dev_vec
        elas_dev = jnp.hstack((elas_dev, jnp.sqrt(3.0) * jnp.log(self.a_vol)))
        
        return self.thermo_elas_class.eval(elas_dev, self.pressure_eos, self.eVref)
                
    def elas_strain_to_cauchy_stress(self, elas_dev_vec):
        kirchoff = self.elas_strain_to_kirchoff_stress(elas_dev_vec)
        return self.inv_det_vol * kirchoff

    def get_elas_strain_state(self, x):
        # Get out elastic strain delta value from solution vector
        elas_delta_dev_vec = jnp.asarray(x[self.ind_e_beg:self.ind_e_end]) * self.elas_scale
        elas_dev_vec_n1 = elas_delta_dev_vec + self.elas_dev_vec_n
        elas_dt_dev_vec = elas_delta_dev_vec * self.inv_delta_time

        return (elas_delta_dev_vec, elas_dev_vec_n1, elas_dt_dev_vec)

    def get_rotation_state(self, x):
        #get out the omega tensor delta value from solution vector
        delta_omega = jnp.asarray(x[self.ind_r_beg:self.ind_r_end]) * self.rot_scale
        crystal_quat_delta = jeu.exp_map_to_quat(delta_omega)
        crystal_quat_n1 = jeu.update_quat_rot(crystal_quat_delta, self.crystal_quat_n)

        crystal_rmat = jeu.quat_to_rmat(crystal_quat_n1)
        crystal_rot_mat5 = jeu.rot_mat_to_rot_mat5(crystal_rmat)

        return (delta_omega, crystal_rmat, crystal_rot_mat5)

    def calc_slip_system_terms(self, crystal_kirchoff_dev):
        # Calculate quantities related to slip system
        # Note not all systems will actually use chia so it might just be a zeros vector
        # We're just combining things here to make our lives a bit less complicated at the
        # cost of efficiency
        chi, schmid_system_p_vecs, schmid_system_q_vecs = self.slip_geom_class.get_PQ_chia(crystal_kirchoff_dev, True)

        # Calculate our resolved shear stress and then slip rates
        rss = self.slip_geom_class.evaluate_RSS(crystal_kirchoff_dev, schmid_system_p_vecs)
        # Eventually we should be able to have the derivative terms calculated for us through AD but for now that's not important
        if self.slip_geom_class.dynamic >= 1:
            rsschi = jnp.hstack([rss, chi])
            slip_rates = self.slip_kinetics_class.eval_slip_rates(rsschi, self.kinetic_vals)
        else:
            slip_rates = self.slip_kinetics_class.eval_slip_rates(rss, self.kinetic_vals)
        # Calculate the plastic slip rate symmetric and skew tensor values
        plastic_def_rate_dev_vecs = jnp.dot(schmid_system_p_vecs, slip_rates)
        plastic_spin_dev_vecs = jnp.dot(schmid_system_q_vecs, slip_rates)

        return (rss, slip_rates, plastic_def_rate_dev_vecs, plastic_spin_dev_vecs)
    
    def get_residual(self, x, args=()):
        # Calculate related elastic strain and lattice rotation terms
        elas_delta_dev_vec, elas_dev_vec_n1, elas_dt_dev_vec = self.get_elas_strain_state(x)
        delta_omega, crystal_rmat, crystal_rot_mat5 = self.get_rotation_state(x)

        # Rotate sample deformation rate tensor and spin vec back to crystal
        crystal_def_dev_vec = jnp.dot(crystal_rot_mat5.T, self.def_dev_vec_samp)
        crystal_spin_dev_vec = jnp.dot(crystal_rmat.T, self.spin_dev_vec_samp)

        # Calculate the deviatoric Kirchoff stress tensor
        crystal_kirchoff_dev = self.elas_strain_to_kirchoff_stress(elas_dev_vec_n1)
        # Calculate quantities related to slip system
        rss, _, plastic_def_rate_dev_vecs, plastic_spin_dev_vecs = self.calc_slip_system_terms(crystal_kirchoff_dev)

        # Can now start calculating other terms related to the residual
        # For the terms related to the change in the change in the omega aka Rmat_dot Rmat term
        # term of our residual
        # Want the edot e - e edot term (a skew matrix) in as a 3x1 value
        # First get out the transformation matrix
        elas_dev_oper_skw = jeu.mat35_da_A_oper_b_d(elas_dev_vec_n1)
        # Need to double check this is what we expect
        elas_dot_elas_dev_skw = jnp.dot(elas_dev_oper_skw, elas_dt_dev_vec) 
        elas_dot_elas_factor = 0.5 * self.inv_a_vol * self.inv_a_vol

        residual_elas = self.epsdot_scale_inv * (self.inv_a_vol * elas_dt_dev_vec + plastic_def_rate_dev_vecs - crystal_def_dev_vec)

        residual_omega = self.rotation_incr_scale_inv * (delta_omega - self.delta_time * (crystal_spin_dev_vec - plastic_spin_dev_vecs + elas_dot_elas_factor * elas_dot_elas_dev_skw))

        return jnp.hstack((residual_elas, residual_omega))

    def get_jacobian(self, x, args=()):
        return jax.jacfwd(self.get_residual, 0)(x)

    def compute_resid_jacobian(self, x, args=()):
        residual = self.get_residual(x)
        jacobian = self.get_jacobian(x)
        return (residual, jacobian)

    def calculate_other_terms(self, x):
        _, elas_dev_vec_n1, _ = self.get_elas_strain_state(x)
        # Calculate the deviatoric Kirchoff stress tensor
        crystal_kirchoff_dev = self.elas_strain_to_kirchoff_stress(elas_dev_vec_n1)
        # Calculate slip system related terms
        rss, slip_rates, plastic_def_rate_dev_vecs, _ = self.calc_slip_system_terms(crystal_kirchoff_dev)

        # Additional factors that we don't really need but could be useful for outside use
        self.slip_rates = jnp.copy(slip_rates)
        self.plastic_disipation_rate_contribution = self.inv_a_vol * jnp.sum(rss * slip_rates)
        self.shear_rate_effective_contribution = jeu.vec_dev_effective(plastic_def_rate_dev_vecs)

def get_response(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class,
                 delta_time, solver_tolerance, def_rate_samp, spin_vec_samp,
                 vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
                 temp_k):

    hist_class = jec.HistClass(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class)

    dmean = -1.0 / 3.0 * (def_rate_samp[0] + def_rate_samp[1] + def_rate_samp[2])

    def_rate_vec7_samp = jnp.asarray([def_rate_samp[0] + dmean,
                                      def_rate_samp[1] + dmean,
                                      def_rate_samp[2] + dmean,
                                      def_rate_samp[3],
                                      def_rate_samp[4],
                                      def_rate_samp[5],
                                      -3.0 * dmean])

    def_dev_vec_samp = jeu.sym_vec_to_vec_dev(def_rate_vec7_samp)

    hard_state_n = hist_class.get_hard_state(history_vec)
    slip_rate_n = hist_class.get_slip_rate(history_vec)
    elas_dev_vec_n = hist_class.get_elas_dev(history_vec)
    # This also normalizes the quats just in-case they weren't ahead of time
    crystal_quat_n = hist_class.get_quats(history_vec)

    # Calculate deviatoric strain energy contribution using trapizodal rule
    half_vol_mid_dt= 0.25 * (vol_ratio_vec[0] + vol_ratio_vec[1]) * delta_time
    beg_dev_strain_energy = half_vol_mid_dt * jeu.inner_prod_sym_vec(stress_vec_pressure, def_rate_vec7_samp)

    # EOS calculations here now
    energy_old = internal_energy[0]
    pressure_old = stress_vec_pressure[-1]

    _, temp_k = eos_class.eval_pressure_temp(vol_ratio_vec[0], energy_old)

    temp_k_new, press_eos, energy_new, bulk_mod_new = jeos.update_simple(eos_class, vol_ratio_vec[1], vol_ratio_vec[3], energy_old, pressure_old)

    def calc_kirchoff_stress():
        a_vol = jnp.power(vol_ratio_vec[1], 1.0 / 3.0)
        inv_a_vol = 1.0 / a_vol
        elas_dev = inv_a_vol * elas_dev_vec_n
        elas_dev = jnp.hstack((elas_dev, jnp.sqrt(3.0) * jnp.log(a_vol))) 
        return thermo_elas_class.eval(elas_dev, press_eos, energy_new)

    crystal_kirchoff_dev = calc_kirchoff_stress()

    # Hardening update using beg of time step values
    chia, _, _ = slip_geom_class.get_PQ_chia(crystal_kirchoff_dev, setvals=True)
    # ignore the nfev value returned as we don't save it anywhere
    _, hard_state_n1 = slip_kinetics_class.update_hardness(hard_state_n, chia, slip_rate_n, delta_time, temp_k)

    # Elastic and lattice rotation updates

    evptn_class = evptnClass(slip_geom_class, slip_kinetics_class, thermo_elas_class,
                             delta_time, vol_ratio_vec[1], energy_new, press_eos,
                             temp_k, hard_state_n1, elas_dev_vec_n, crystal_quat_n,
                             def_dev_vec_samp, spin_vec_samp)

    x0 = jnp.zeros(jec.NWVEC + jec.NTVEC)
    # res = root(evptn_class.compute_resid_jacobian, x0, jac=True, method='hybr', tol=1e-8)
    # solver = snls.SNLSTrDlDenseG(evptn_class.compute_resid_jacobian, xtolerance=solver_tolerance, ndim=x0.shape[0])
    # solver.delta_control.deltaInit = 1.0
    # status, xs = solver.solve()
    # xs = res.x
    solver = optx.Dogleg(rtol=1e-6, atol=1e-8)
    sol = optx.root_find(evptn_class.get_residual, solver=solver, y0=x0, args=())
    xs = sol.value
    evptn_class.calculate_other_terms(xs)

    elas_dev_vec_n1, crystal_quat_n1 = evptn_class.get_state_from_x(xs)
    slip_rate_n1 = jnp.copy(evptn_class.slip_rates)

    shear_eff = hist_class.get_shear_eff(history_vec)
    shear_rate_eff = evptn_class.shear_rate_effective_contribution
    shear_eff += shear_rate_eff * delta_time

    def_effective = jeu.vec_dev_effective(def_dev_vec_samp)

    flow_strength = jax.lax.cond(
        def_effective > jec.DBL_TINY_SQRT,
        lambda: evptn_class.plastic_disipation_rate_contribution / def_effective,
        lambda: evptn_class.hard_scale
    )

    solver_iters = sol.stats["num_steps"]

    cauchy_crystal = evptn_class.elas_strain_to_cauchy_stress(elas_dev_vec_n1)

    rmat_n1 = jeu.quat_to_rmat(crystal_quat_n1)
    rmat_m5 = jeu.rot_mat_to_rot_mat5(rmat_n1)

    cauchy_samp_dev_vec = jnp.dot(rmat_m5, cauchy_crystal[0:5])
    cauchy_samp = jnp.hstack((cauchy_samp_dev_vec, cauchy_crystal[-1]))
    stress_vec_pressure_n1 = jeu.dev_vec_to_sym_vec(cauchy_samp)

    dev_strain_energy = beg_dev_strain_energy + half_vol_mid_dt * jeu.inner_prod_sym_vec(stress_vec_pressure_n1, def_rate_vec7_samp)

    sdd = jnp.asarray([bulk_mod_new, thermo_elas_class.shear_mod])

    energy_new += dev_strain_energy

    internal_energy_n1 = jnp.asarray([energy_new])

    crystal_quat_n1 = jax.lax.cond(
        jnp.sum(crystal_quat_n * crystal_quat_n1) < 0.0,
        lambda: crystal_quat_n1 * -1.0,
        lambda: crystal_quat_n1
    )

    history_update = hist_class.pack_history_vars(elas_dev_vec_n1, crystal_quat_n1, hard_state_n1, slip_rate_n1, shear_rate_eff, shear_eff, flow_strength, solver_iters)

    stress_vec = stress_vec_pressure_n1[0:-1]
    stress_vec = stress_vec.at[0:3].set(stress_vec[0:3] - stress_vec_pressure_n1[-1])

    return (stress_vec, (history_update, internal_energy_n1, temp_k, sdd, jnp.copy(xs)))

# Due to the above nonlinear solver consistently giving NANs during the backpropagation step to get
# out the material tangent stiffness matrix, we needed to resort to essentially duplicating portions
# of the above so that we could get out correct values of the material tangent stiffness matrix when doing AD calcs :/
def get_response_mtan(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class,
                      delta_time, def_rate_samp, spin_vec_samp,
                      vol_ratio_vec, history_vec, ie_peos_tk_hard_tup, xs):

    hist_class = jec.HistClass(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class)

    dmean = -1.0 / 3.0 * (def_rate_samp[0] + def_rate_samp[1] + def_rate_samp[2])

    def_rate_vec7_samp = jnp.asarray([def_rate_samp[0] + dmean,
                                      def_rate_samp[1] + dmean,
                                      def_rate_samp[2] + dmean,
                                      def_rate_samp[3],
                                      def_rate_samp[4],
                                      def_rate_samp[5],
                                      -3.0 * dmean])

    def_dev_vec_samp = jeu.sym_vec_to_vec_dev(def_rate_vec7_samp)

    slip_rate_n = hist_class.get_slip_rate(history_vec)
    elas_dev_vec_n = hist_class.get_elas_dev(history_vec)
    # This also normalizes the quats just in-case they weren't ahead of time
    crystal_quat_n = hist_class.get_quats(history_vec)

    energy_new, press_eos, temp_k, hard_state_n1 = ie_peos_tk_hard_tup

    def calc_kirchoff_stress():
        a_vol = jnp.power(vol_ratio_vec[1], 1.0 / 3.0)
        inv_a_vol = 1.0 / a_vol
        elas_dev = inv_a_vol * elas_dev_vec_n
        elas_dev = jnp.hstack((elas_dev, jnp.sqrt(3.0) * jnp.log(a_vol))) 
        return thermo_elas_class.eval(elas_dev, press_eos, energy_new)

    crystal_kirchoff_dev = calc_kirchoff_stress()

    # Hardening update using beg of time step values
    chia, _, _ = slip_geom_class.get_PQ_chia(crystal_kirchoff_dev, setvals=True)
    # Elastic and lattice rotation updates

    evptn_class = evptnClass(slip_geom_class, slip_kinetics_class, thermo_elas_class,
                             delta_time, vol_ratio_vec[1], energy_new, press_eos,
                             temp_k, hard_state_n1, elas_dev_vec_n, crystal_quat_n,
                             def_dev_vec_samp, spin_vec_samp)

    # Technically not correct but this does allow JAX to get out more or less the correct
    # material tangent stiffness matrix. If we could tell JAX to only worry about the last
    # iteration of our various nonlinear solver updates for the jacobian calcs then this
    # wouldn't be an issue...
    resid, jacob = evptn_class.compute_resid_jacobian(xs)
    xsol = xs - jnp.linalg.solve(jacob, resid)

    resid = evptn_class.get_residual(xsol)

    evptn_class.calculate_other_terms(xsol)

    elas_dev_vec_n1, crystal_quat_n1 = evptn_class.get_state_from_x(xsol)

    cauchy_crystal = evptn_class.elas_strain_to_cauchy_stress(elas_dev_vec_n1)

    rmat_n1 = jeu.quat_to_rmat(crystal_quat_n1)
    rmat_m5 = jeu.rot_mat_to_rot_mat5(rmat_n1)

    cauchy_samp_dev_vec = jnp.dot(rmat_m5, cauchy_crystal[0:5])
    cauchy_samp = jnp.hstack((cauchy_samp_dev_vec, cauchy_crystal[-1]))
    stress_vec_pressure_n1 = jeu.dev_vec_to_sym_vec(cauchy_samp)

    stress_vec = stress_vec_pressure_n1[0:-1]
    stress_vec = stress_vec.at[0:3].set(stress_vec[0:3] - stress_vec_pressure_n1[-1])

    return stress_vec

if __name__ == "__main__":

    params = {}
    #general parameters
    params["sol_tolerance"] = 1e-10
    # thermo elast n parameters
    params["C11"] = 300e-2
    params["C12"] = 100e-2
    params["C44"] = 100e-2
    #EOS parameters minus bulk modulus

    params["init_density"] = 3.0
    params["cvav"] = 2.0e-5
    params["eos_gamma"] = 1.7
    params["eos_cold_energy_0"] = -params["cvav"] * 300.0
    # dummy parameter for now will update after thermo elast n class created
    # and then use their bulk modulus calculation
    params["bulk_modulus_0"] = -1.0
    #slip geometry parameters
    params["num_slip_systems"] = 12
    # slip kinetics and hardening parameters
    params["slip_kin_nonlinear"] = False
    params["shear_mod"] = 1.0 
    params["slip_kin_exp_m"] = 0.01
    params["slip_kin_gamma_0_w"] = 1.0
    params["slip_kin_h0"] = 200e-5 
    params["slip_kin_crss0"] = 100e-5
    params["slip_kin_crss_sat"] = 400e-5
    params["slip_kin_voce_exp_n"] = 1.0
    params["slip_kin_voce_exp_m_sat"] = 0.05
    params["slip_kin_voce_gamma_sat_0"] = 1.0e-6

    thermo_elas_class = jtelas.thermoElastCubic(params)
    # update our bulk modulus value based on what thelc calculated
    params["bulk_modulus_0"] = thermo_elas_class.bulk_mod
    # Isothermal case
    eos_class = jeos.eosSimple(params, False) 
    slip_geom_class = jslgeo.SlipGeomFCC(params)
    slip_kinetics_class = jslkin.SlipKineticVocePowerLaw(params)

    hist_class = jec.HistClass(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class)
    num_hist = hist_class.num_hist

    # Need to come up with a better way to do this later on to set-up history variables...
    elas_dev = jnp.zeros(jec.NTVEC)
    quats = jnp.asarray([1.0, 0.0, 0.0, 0.0])
    hard_state = jnp.asarray([params["slip_kin_crss0"]])
    slip_rate  = jnp.zeros(slip_geom_class.num_slip_systems)
    shear_rate_eff = 0.0
    shear_eff = 0.0
    flow_strength = 0.0
    solver_iters = 0.0

    history_vec = hist_class.pack_history_vars(elas_dev, quats, hard_state, slip_rate, shear_rate_eff, 
                                               shear_eff, flow_strength, solver_iters)

    jax.debug.print("history vec {}", history_vec)

    #test conditions
    def_rate_vec7_samp = jnp.asarray([-0.5, -0.5, 1.0, 0.001, 0.001, 0.001]) * jnp.sqrt(2.0/3.0)
    spin_vec_samp = jnp.asarray([0.0, 0.0, 0.5])
    vol_ratio_vec = jnp.asarray([1.0, 1.0, 0.0, 0.0])

    delta_time = 1e-1
    internal_energy = jnp.asarray([0.0])
    stress_vec_pressure = jnp.zeros((jec.NSVP))
    temp_k = 300.

    get_response_jit = jax.jit(get_response, static_argnums=(0, 1, 2, 3, 5))

    stress_vec_pressure_n1, others = get_response_jit(
                 slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class,
                 delta_time, params["sol_tolerance"], def_rate_vec7_samp, spin_vec_samp,
                 vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
                 temp_k
                 )

    history_update, internal_energy_n1, temp_k_n1, sdd = others

    # jacobians, others = jax.jacrev(get_response, argnums=6, has_aux=True)(
    #              slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class,
    #              delta_time, params["sol_tolerance"], def_rate_vec7_samp, spin_vec_samp,
    #              vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
    #              temp_k
    #              )
    # Add back the bulk contribution for the tangent stiffness matrix
    # jacob_bulk = np.zeros((6, 6))
    # jacob_bulk[-1,-1] = 3.0 * sdd[0]
    # jacob_bulk *= delta_time
    # jacob_bulk = jeu.mtan_conv_sd_svec(jacob_bulk, True)

    # jacob_np = np.asarray(jacobians) + jacob_bulk
    # jacob_np[3:-1, 3:-1] *= 0.5

    print(stress_vec_pressure_n1)
    print()
    print(internal_energy_n1)
    print()
    print(temp_k)
    print()
    print(sdd)
    print("Slip Rates")    
    jax.debug.print("{}", hist_class.get_slip_rate(history_update))
    print("Deviatoric crystal elastic strain")
    jax.debug.print("{}", hist_class.get_elas_dev(history_update))
    print("Lattice quaternions")
    jax.debug.print("{}", hist_class.get_quats(history_update))
    print("Number of function evaluations")
    jax.debug.print("{}", history_update[hist_class.ind_hist_num_func_evals])
