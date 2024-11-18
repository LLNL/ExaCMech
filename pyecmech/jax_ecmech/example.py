#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 08:22:12 2023

@author: carson16
"""

import numpy as np
import jax

#Can't save the jax compilation between runs yet :/
# from jax.experimental.compilation_cache import compilation_cache as cc
# cc.set_cache_dir("./jax-cache")

import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

import jax_evptn_wrap as jevptnw
import jax_ecmech_const as jecm

class JaxProb:
    def __init__(self, params):
        self.prob = jevptnw.evptnWrapClass(params)

    def init_history_vec(self):
        return self.prob.init_history_vec()

    def get_history_info(self):
        return self.prob.get_history_info()

    def solve(self, delta_time, def_rate_samp, spin_vec_samp, vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec, temp_k, need_mtan=False):
        '''
            Copied from example.py file for their solve case and need to update for batch solves
            as currently this only deals with single point solves...

            Solve does a per time step solve of the material update for all points inputted.
            A few things to note:
            def_rate_samp has dimensions npts x self.nsvec input
            spin_vec_samp has dimesnions npts x self.nwvec input
            vol_ratio_vec has dimensions npts x self.nvr input
            internal_energy has dimensions npts x self.ne input/output
            stress_vec_pressure has dimensions npts x self.nsvp input/output
            history_vec has dimensions npts x self.nhist input/output
            temp_k has dimensions npts x 1 input/output
            need_mtan - Set to True if you need the tangent stiffness matrix

            returns a tuple of:
            stress_vec_pressure_n1, history_update, internal_energy_n1, temp_k, sdd, and potentially mtan

            If you pass in 1D arrays we will promote them to 2D arrays.
        '''
        # Need to update the underlying models to get rid of all the if statements but eventually we should
        # just be able to use a jax vectorized / jit'd version of things...

        bdt = jnp.atleast_1d(delta_time)
        drs = jnp.atleast_2d(def_rate_samp)
        svs = jnp.atleast_2d(spin_vec_samp)
        vrv = jnp.atleast_2d(vol_ratio_vec)
        ie  = jnp.atleast_2d(internal_energy)
        svp = jnp.atleast_2d(stress_vec_pressure)
        hiv = jnp.atleast_2d(history_vec)
        tk  = jnp.atleast_1d(temp_k)

        npts = tk.shape[0]

        stress_vec_pressure_n1 = np.zeros_like(svp)
        history_update = jnp.zeros_like(hiv)
        sdd = np.zeros((npts, jecm.NSDD))
        internal_energy_n1 = np.zeros_like(ie)
        temp_k_n1 = np.zeros_like(temp_k)

        if need_mtan:
            jacob_np = np.zeros((npts, 6, 6))
        else:
            jacob_np = None
        
        for ipts in range(npts):
            svp_n1, hu, ie_n1, temp, sdd_t, targ = self.prob.solve(delta_time, drs[ipts, :], svs[ipts, :], vrv[ipts, :], ie[ipts, :], svp[ipts, :], hiv[ipts, :], tk[ipts], need_mtan)

            stress_vec_pressure_n1[ipts, :] = np.asarray(svp_n1)
            # I still don't understand why some of these values I can easily convert to numpy arrays and others can't :(
            history_update = history_update.at[ipts, :].set(hu)
            temp_k_n1[ipts] = np.asarray(temp)
            sdd[ipts, :] = np.asarray(sdd_t)

            if need_mtan:
                jacob_np[ipts, :, :] = np.asarray(targ)

        return (stress_vec_pressure_n1, history_update, internal_energy_n1, temp_k, sdd, jacob_np)


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

    params["slip_system_geometry"] = "fcc"
    params["thermo_elast_class"] = "thermoElastCubic"
    params["eos_class_isothermal"] = "false"
    params["slip_kinetics_hardening_class"] = "voce_pl"

    prob = JaxProb(params)

    history_vec = prob.init_history_vec()

    jax.debug.print("history vec {}", history_vec)

    npts = 1
    nsteps = 41
    #test conditions
    # Note the deformation rate here is the full tensor but as a 6d vector
    # Internally, we will decompose things into a deviatoric and volumetric component
    def_rate_samp = jnp.asarray([-0.5, -0.5, 1.0, 0.001, 0.001, 0.001]) * jnp.sqrt(2.0/3.0)
    spin_vec_samp = jnp.asarray([0.0, 0.0, 0.5])
    vol_ratio_vec = jnp.asarray([1.0, 1.0, 0.0, 0.0])

    delta_time = 1e-1
    internal_energy = jnp.asarray([0.0])
    stress_vec_pressure = jnp.zeros((jecm.NSVP))
    temp_k = np.ones(1) * 300.

    # Making data 2d here to test batch solves could have also just called
    # jnp.atleast_2d(arr) here instead
    def_rate_samp = jnp.tile(def_rate_samp, (npts, 1))
    spin_vec_samp = jnp.tile(spin_vec_samp, (npts, 1))
    internal_energy = jnp.tile(internal_energy, (npts, 1))
    stress_vec_pressure = jnp.tile(stress_vec_pressure, (npts, 1))
    history_vec = jnp.tile(history_vec, (npts, 1))
    vol_ratio_vec = jnp.tile(vol_ratio_vec, (npts, 1))
    temp_k = jnp.repeat(temp_k, npts)

    # Note the last value returned here is the material tangent stiffness matrix
    # but we're just going to ignore that here
    for i in range(nsteps):
        dmean = def_rate_samp[:, 0] + def_rate_samp[:, 1] + def_rate_samp[:, 2]

        vol_ratio_vec = vol_ratio_vec.at[:, 0].set(vol_ratio_vec[:, 1])
        vol_ratio_vec = vol_ratio_vec.at[:, 1].set(vol_ratio_vec[:, 0] * jnp.exp(dmean * delta_time))
        vol_ratio_vec = vol_ratio_vec.at[:, 3].set(vol_ratio_vec[:, 1] - vol_ratio_vec[:, 0])
        vol_ratio_vec = vol_ratio_vec.at[:, 2].set(vol_ratio_vec[:, 3] / (delta_time * 0.5 * (vol_ratio_vec[:, 0] + vol_ratio_vec[:, 1])))

        stress_vec_pressure, history_vec, internal_energy, temp_k, sdd, _ = prob.solve(
                    delta_time, def_rate_samp, spin_vec_samp, vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec, temp_k)

        print(stress_vec_pressure[:,0])

    print("Deviatoric Stress + pressure:")
    print(stress_vec_pressure)
    print("Internal energy")
    print(internal_energy)
    print("Temperature (K)")
    print(temp_k)
    print("SDD array")
    print(sdd)
    print("Slip Rates")
    mean_history_update = jnp.mean(history_vec, axis=0) 
    hc = prob.prob.hist_class
    jax.debug.print("{}", hc.get_slip_rate(mean_history_update))
    print("Deviatoric crystal elastic strain")
    jax.debug.print("{}", hc.get_elas_dev(mean_history_update))
    print("Lattice quaternions")
    jax.debug.print("{}", hc.get_quats(mean_history_update))
    print("Number of function evaluations")
    jax.debug.print("{}", mean_history_update[hc.ind_hist_num_func_evals])




