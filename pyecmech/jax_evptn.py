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

import jax_ecmech_util as jeu
import jax_ecmech_const as jec


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
        self.a_vol = jnp.pow(det_vol, 1.0 / 3.0)
        self.inv_a_vol = 1.0 / self.a_vol
        
        self.hard_scale, self.kinetic_vals = self.slip_kinetics_class.get_vals(self.pressure_eos, self.temp_k, self.hard_state)

        adot_ref = self.slip_kinetics_class.get_fixed_ref_rate(self.kinetic_vals)
        eff = jnp.norm(def_dev_vec_samp)
        
        self.epsdot_scale_inv = jax.lax.cond(
            eff < (jnp.sqrt(jnp.finfo(jnp.float64).eps) * adot_ref),
            lambda: 1.0 / adot_ref,
            lambda: jnp.min(1.0 / eff, 1.0e6 * delta_time)
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
        elas_dev = jnp.c_[elas_dev, jnp.sqrt(3.0) * jnp.log(self.a_vol)]
        
        return self.thermo_elas_class.eval(elas_dev, self.temp_k, self.pressure_eos, self.eVref)
                
    def elas_strain_to_cauchy_stress(self, elas_dev_vec):
        kirchoff = self.elas_strain_to_kirchoff_stress(elas_dev_vec)
        return self.inv_det_vol * kirchoff
    
    def get_residual(self, x):
        # Get out elastic strain delta value from solution vector
        elas_delta_dev_vec = x[self.ind_e_beg:self.ind_e_end] * self.elas_scale
        elas_dev_vec_n1 = elas_delta_dev_vec + self.elas_dev_vec_n
        elas_dt_dev_vec = elas_delta_dev_vec * self.inv_delta_time

        #get out the omega tensor delta value from solution vector
        delta_omega = x[self.ind_r_beg:self.ind_r_end] * self.r_scale
        crystal_quat_delta = jeu.exp_map_to_quat(omega_spin)
        crystal_quat_n1 = jeu.update_quat_rot(crystal_quat_delta, self.crystal_quat_n)

        crystal_rmat = jeu.quat_to_rmat(crystal_quat_n1)
        crystal_rot_mat5 = jeu.rot_mat_to_rot_mat5(crystal_rmat)

        # Rotate sample deformation rate tensor and spin vec back to crystal
        crystal_def_dev_vec = jnp.dot(crystal_rot_mat5.T, self.def_dev_vec_samp)
        crystal_spin_dev_vec = jnp.dot(crystal_rmat.T, self.spin_dev_vec_samp)

        # Calculate the deviatoric Kirchoff stress tensor
        crystal_kirchoff_dev = self.elas_strain_to_kirchoff_stress(elas_dev_vec_n1)

        # Calculate quantities related to slip system
        # Note not all systems will actually use chia so it might just be a zeros vector
        # We're just combining things here to make our lives a bit less complicated at the
        # cost of efficiency
        schmid_system_p_vecs, schmid_system_q_vecs, chia = self.slip_geom_class.get_PQ_chia(crystal_kirchoff_dev)

        # Calculate our resolved shear stress and then slip rates
        rss = self.slip_geom_class.evaluate_RSS(crystal_kirchoff_dev)
        # Eventually we should be able to have the derivative terms calculated for us through AD but for now that's not important
        self.slip_rates = self.slip_kinetics_class.eval_slip_rates(rss, self.kinetic_vals)
        # Calculate the plastic slip rate symmetric and skew tensor values
        plastic_def_rate_dev_vecs = jnp.dot(schmid_system_p_vecs, self.slip_rates)
        plastic_spin_dev_vecs = jnp.dot(schmid_system_q_vecs, self.slip_rates)

        # Additional factors that we don't really need but could be useful for outside use
        self.plastic_disipation_rate_contribution = self.inv_a_vol * jnp.sum(rss * self.slip_rates)
        self.shear_rate_effective_contribution = jeu.vec_dev_effective(plastic_def_rate_dev_vecs)

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

        return jnp.c_[residual_elas, residual_omega]

    def get_jacobian(self, x):
        return jax.jacfwd(self.get_residual, argnums=1)(x)

    def compute_resid_jacobian(self, x):
        residual = self.get_residual(x)
        jacobian = self.get_jacobian(x)

        return (residual, jacobian)


def get_response(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class,
                 delta_time, solver_tolerance, def_rate_vec7_samp, spin_vec_samp,
                 vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
                 temp_k):

    hist_class = jec.HistClass(slip_geom_class, slip_kinetics_class, thermo_elas_class, eos_class)

    def_dev_vec_samp = jeu.sym_vec_to_vec_dev(def_rate_vec7_samp)

    hard_state_n = hist_class.get_hard_state(hist)
    slip_rate_n = hist_class.get_slip_rate(hist)
    elas_dev_vec_n = hist_class.get_elas_dev(hist)
    # This also normalizes the quats just in-case they weren't ahead of time
    crystal_quat_n = hist_class.get_quats(hist)

    # Calculate deviatoric strain energy contribution using trapizodal rule
    half_vol_mid_dt= 0.25 * (vol_ratio_vec[0] + vol_ratio_vec[1]) * delta_time
    beg_dev_strain_energy = half_vol_mid_dt * jeu.inner_prod_sym_vec(stress_vec_pressure, def_rate_vec7_samp)

    # EOS calculations here now
    energy_old = internal_energy[0]
    pressure_old = stress_vec_pressure[-1]

    temp_k, junk = eos_class.eval_pressure_temp(vol_ratio_vec[0], energy_old)

    temp_k_new, press_eos, energy_new, bulk_mod_new, junk1, junk2, junk3 = jeeos.update_simple(eos_class, vol_ratio_vec[1], vol_ratio_vec[3], energy_old, pressure_old)

    # Hardening update using beg of time step values
    schmid_system_p_vecs, schmid_system_q_vecs, chia = self.slip_geom_class.get_PQ_chia(crystal_kirchoff_dev)

    hard_state_n1 = slip_kinetics_class.update_hardness(hard_state_n, delta_time, slip_rate_n, chia, temp_k)

    # Elastic and lattice rotation updates

    evptn_class = evptnClass(slip_geom_class, slip_kinetics_class, thermo_elas_class,
                             delta_time, vol_ratio_vec[1], energy_new, press_eos,
                             temp_k, hard_state_n1, elas_dev_vec_n, crystal_quat_n,
                             def_dev_vec_samp, spin_dev_vec_samp)

    x0 = jnp.zeros(jec.NWVEC + jec.NTVEC)
    res = sciop.root(evptn_class.compute_resid_jacobian, x0, jac=True, method='hybr', tol=1e-8)

    elas_dev_vec_n1, crystal_quat_n1 = evptn_class.get_state_from_x(res.x)
    slip_rate_n1 = jnp.copy(evptn_class.slip_rates)

    shear_eff = hist_class.get_shear_eff(hist)
    shear_rate_eff = evptn_class.shear_rate_effective_contribution
    shear_eff += shear_rate_eff * delta_time

    def_effective = jeu.vec_dev_effective(def_dev_vec_samp)

    if def_effective > jec.DBL_TINY_SQRT:
        flow_strength = evptn_class.plastic_disipation_rate_contribution / def_effective
    else:
        flow_strength = self.hard_scale

    solver_iters = res.nit

    cauchy_crystal = evptn_class.elas_strain_to_cauchy_stress(elas_dev_vec_n1)

    rmat_n1 = jeu.quat_to_rmat(crystal_quat_n1)
    rmat_m5 = jeu.rot_mat_to_rot_mat5(rmat_n1)

    cauchy_samp_dev_vec = jnp.dot(rmat_m5, cauchy_crystal[0:6])
    cauchy_samp = jnp.c_[cauchy_crystal_dev_vec, cauchy_crystal[-1]]
    stress_vec_pressure_n1 = dev_vec_to_sym_vec(cauchy_samp)

    dev_strain_energy = beg_dev_strain_energy + half_vol_mid_dt * jeu.inner_prod_sym_vec(stress_vec_pressure_n1, def_rate_vec7_samp)

    sdd = jnp.asarray([bulk_mod_new, thermo_elas_class.get_shear_mod(temp_k, press_eos, energy_new)])

    energy_new += dev_strain_energy

    internal_energy_n1 = jnp.asarray([energy_new])

    if jnp.sum(crystal_quat_n * crystal_quat_n1) < 0.0:
        crystal_quat_n1 *= -1.0

    history_update = hist_class.pack_history_vars(elas_dev_vec_n1, crystal_quat_n1, hard_state_n1, slip_rate_n1, shear_rate_eff, shear_eff, flow_strength, solver_iters)

    return (stress_vec_pressure_n1, history_update, internal_energy_n1, temp_k, sdd)






