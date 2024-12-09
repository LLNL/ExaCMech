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

import jax_ecmech_util as jeu
import jax_ecmech_const as jec

import jax_eos as jeos
import jax_slip_geom as jslgeo
import jax_slip_kinetics as jslkin
import jax_thermo_elastn as jtelas
import jax_evptn as jevptn

np.set_printoptions(linewidth=np.inf)

class evptnWrapClass:
    def __init__(
                 self,
                 params
                 ):

        match params["slip_system_geometry"]:
            case "fcc" | "FCC":
                self.crystal_symmetry = "cubic"
                self.slip_geom_class = jslgeo.SlipGeomFCC(params)
            case "bcc" | "BCC":
                self.crystal_symmetry = "cubic"
                self.slip_geom_class = jslgeo.SlipGeomBCC(params)
            case "bcc_pencil":
                self.crystal_symmetry = "cubic"
                self.slip_geom_class = jslgeo.SlipGeomBCCPencil(params)
            case _ :
                val = params["slip_system_geometry"]
                raise ValueError(f"A slip_system_geometry value was not provided {val}")

        match params["thermo_elast_class"]:
            case "ThermoElasicCubic" | "thermoElastCubic":
                if self.crystal_symmetry != "cubic":
                    val_ssg = params["slip_system_geometry"]
                    val_tec = params["thermo_elast_classs"]
                    raise ValueError(f"Trying to use incompatible crystal symmetry: {self.crystal_symmetry} as determined by slip_system_geometry: {val_ssg}  and thermo_elast_class: {val_tec}")
                self.thermo_elas_class = jtelas.thermoElastCubic(params)
            case _ :
                val = params["thermo_elast_class"]
                raise ValueError(f"A thermo_elast_class value was not provided {val}")

        # update our bulk modulus value based on what thermo_elas_class calculated
        params["bulk_modulus_0"] = self.thermo_elas_class.bulk_mod
        if not "shear_mod" in params:
                    params["shear_mod"] = self.thermo_elas_class.shear_mod

        match params["eos_class_isothermal"]:
            case "isothermal" | "true" | True:
                # Isothermal case
                self.eos_class = jeos.eosSimple(params, True)
            case _ :
                self.eos_class = jeos.eosSimple(params, False)


        match params["slip_kinetics_hardening_class"]:
            case "voce_pl" | "Voce_PL":
                self.slip_kinetics_class = jslkin.SlipKineticVocePowerLaw(params)
            case "oro_dd":
                self.slip_kinetics_class = jslkin.SlipKineticOrowanD(params)
            case "km_bal_dd":
                self.slip_kinetics_class = jslkin.SlipKineticMTSKocksMecking(params)
            case "bcc_md":
                self.slip_kinetics_class = jslkin.SlipKineticBCCMD(params)
            case _ :
                val = params["slip_kinetics_hardening_class"]
                raise ValueError(f"A slip_kinetics_hardening_class value was not provided {val}")

        self.solver_tolerance = params["sol_tolerance"]

        self.hist_class = jec.HistClass(self.slip_geom_class, self.slip_kinetics_class, self.thermo_elas_class, self.eos_class)
        self.num_hist = self.hist_class.num_hist

        # Plain JAX jit funcs
        self.get_response_jit = jax.jit(self.get_response_jittable)
        self.mtan_jit = jax.jit(jax.jacrev(self.get_response_mtan_jittable, argnums=1, has_aux=False))
        # Still a WIP to get all the necessary things ported to JAX idioms so that
        # we can have vectorized calls
        # So this does at least appear to work as the code doesn't crash...
        # No idea if it actually works though...
        # self.batch_solve = jax.vmap(self.solve)

    def get_response_jittable(self, 
                 delta_time, def_rate_samp, spin_vec_samp,
                 vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
                 temp_k):
        return jevptn.get_response(self.slip_geom_class, self.slip_kinetics_class, self.thermo_elas_class, self.eos_class, 
                 delta_time, self.solver_tolerance, def_rate_samp, spin_vec_samp,
                 vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
                 temp_k)

    def get_response_mtan_jittable(self, 
                 delta_time, def_rate_samp, spin_vec_samp,
                 vol_ratio_vec, history_vec, ie_peos_tk_hard_tup, sols):
        return jevptn.get_response_mtan(self.slip_geom_class, self.slip_kinetics_class, self.thermo_elas_class, self.eos_class, 
                 delta_time, def_rate_samp, spin_vec_samp,
                 vol_ratio_vec, history_vec, ie_peos_tk_hard_tup, sols)

    def init_history_vec(self, elas_dev=None, quats=None, hard_state=None, slip_rate=None, shear_rate_eff=None, shear_eff=None, flow_strength=None):
        if elas_dev is None:
            elas_dev = jnp.zeros(jec.NTVEC)
        if quats is None:
            quats = jnp.asarray([1.0, 0.0, 0.0, 0.0])
        if hard_state is None:
            # currently only sub-module which would have non-trivial values here
            _, hard_state, _, _ = self.slip_kinetics_class.get_history_info(list(), list(), list(), list())
            hard_state = jnp.asarray(hard_state)
        if slip_rate is None:
            slip_rate  = jnp.zeros(self.slip_geom_class.num_slip_systems)
        if shear_rate_eff is None:
            shear_rate_eff = 0.0
        if shear_eff is None:
            shear_eff = 0.0
        if flow_strength is None:
            flow_strength = 0.0
        solver_iters = 0.0

        history_vec = self.hist_class.pack_history_vars(
            elas_dev, quats, hard_state, slip_rate, shear_rate_eff, 
            shear_eff, flow_strength, solver_iters)
        
        return history_vec

    def get_history_info(self, names=list(), init=list(), plot=list(), state=list()):

        names.append("eff_plastic_def_rate"); init.append(0.); plot.append(True); state.append(True)
        names.append("equiv_pl_strain"); init.append(0.); plot.append(True); state.append(True)
        names.append("flow_strength"); init.append(0.); plot.append(True); state.append(False)
        names.append("num_func_eval"); init.append(0); plot.append(True); state.append(False)

        for itvec in range(jec.NTVEC):
            name = "xtal_elas_dev_strain_" + str(itvec)
            names.append(name); init.append(0.); plot.append(True); state.append(True)

        names.append(["lattice_quat_0", "lattice_quat_1", "lattice_quat_2", "lattice_quat_3"])
        init.append([1.0, 0.0, 0.0, 0.0])
        plot.append([True, True, True, True])
        state.append([True, True, True, True])

        names, init, plot, state = self.slip_kinetics_class.get_history_info(names, init, plot, state)

        for islip in range(self.slip_geom_class.num_slip_systems):
            name = "shear_rate_" + str(islip)
            names.append(name)
            init.append(0)
            plot.append(True)
            state.append(True)

        return (names, init, plot, state)

    def get_parameters(self):
        params = {}

        params = self.slip_geom_class.get_parameters(params)
        params = self.eos_class.get_parameters(params)
        params = self.thermo_elas_class.get_parameters(params)
        params = self.slip_kinetics_class.get_parameters(params)

        return params

    def mtan_calc(self, args=()):
        delta_time, def_rate_samp, spin_vec_samp, vol_ratio_vec, history_vec, ie_peos_tk_hard_tup, sols, sdd = args
        jacobians = self.mtan_jit(
                delta_time, def_rate_samp, spin_vec_samp,
                vol_ratio_vec, history_vec,
                ie_peos_tk_hard_tup, sols
            )

        jacob_bulk = np.zeros((6, 6))
        jacob_bulk[-1,-1] = 3.0 * sdd[0]
        jacob_bulk *= delta_time
        jacob_bulk = jeu.mtan_conv_sd_svec(jacob_bulk, True)
        jacobians = jnp.where(jnp.abs(jacobians) > 1e-16, jacobians, 0.0)
        jacob_np = np.asarray(jacobians) + jacob_bulk
        jacob_np = np.where(np.abs(jacob_np) > 1e-16, jacob_np, 0.0)
        jacob_np[0:6, 3:6] *= 0.5

        return jacob_np

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
        # In order to make sure we get out the right derivative information later on if needed,
        # we return this as the full Cauchy tensor rather than the 6d deviatoric + pressure variation
        # stress_vec, others = jevptn.get_response(
        #             self.slip_geom_class, self.slip_kinetics_class, self.thermo_elas_class, self.eos_class,
        #             delta_time, self.solver_tolerance, def_rate_samp, spin_vec_samp,
        #             vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
        #             temp_k
        #         )

        stress_vec, others = self.get_response_jit(
                        jnp.array(delta_time), def_rate_samp, spin_vec_samp,
                        vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec,
                        jnp.array(temp_k)
                    )

        pressure = -jnp.sum(stress_vec[0:3]) / 3.0
        stress_vec = stress_vec.at[0:3].set(stress_vec[0:3] + pressure)
        stress_vec_pressure_n1 = jnp.hstack((stress_vec, pressure))
        history_update, internal_energy_n1, temp_k, sdd, sols = others

        jacob_np = None
        if need_mtan:
            ie_peos_tk_hard_tup = (internal_energy_n1, pressure, temp_k, self.hist_class.get_hard_state(history_update))
            args = (delta_time, def_rate_samp, spin_vec_samp, vol_ratio_vec, history_vec, ie_peos_tk_hard_tup, sols, sdd)
            jacob_np = self.mtan_calc(args)

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

    evptn_wc = evptnWrapClass(params)
    quats = np.asarray([1.0, 0.0, 0.0, 0.0])
    history_vec = evptn_wc.init_history_vec(quats=quats)

    jax.debug.print("history vec {}", history_vec)

    #test conditions
    # Note the deformation rate here is the full tensor but as a 6d vector
    # Internally, we will decompose things into a deviatoric and volumetric component
    def_rate_samp = jnp.asarray([-0.5, -0.5, 1.0, 0.001, 0.001, 0.001]) * jnp.sqrt(2.0/3.0)
    # def_rate_vec7_samp = jnp.asarray([-0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 0.0]) * jnp.sqrt(2.0/3.0)
    spin_vec_samp = jnp.asarray([0.0, 0.0, 0.5])
    vol_ratio_vec = jnp.asarray([1.0, 1.0, 0.0, 0.0])

    delta_time = 1e-1
    internal_energy = jnp.asarray([0.0])
    stress_vec_pressure = jnp.zeros((jec.NSVP))
    temp_k = 300.

    # Note the last value returned here is the material tangent stiffness matrix
    # It is only calculated if the user asks us to
    stress_vec_pressure_n1, history_update, internal_energy_n1, temp_k, sdd, mtan = evptn_wc.solve(
                 delta_time, def_rate_samp, spin_vec_samp, vol_ratio_vec, internal_energy, stress_vec_pressure, history_vec, temp_k, True)

    # Still working on getting all the conditionals into a JAX friendly manner so we can vectorize and JIT
    # compile things if need be...
    # bdt = np.atleast_1d(delta_time)
    # drs = np.atleast_2d(def_rate_samp)
    # svs = np.atleast_2d(spin_vec_samp)
    # vrv = np.atleast_2d(vol_ratio_vec)
    # ie  = np.atleast_2d(internal_energy)
    # svp = np.atleast_2d(stress_vec_pressure)
    # hiv = np.atleast_2d(history_vec)
    # tk  = np.atleast_1d(temp_k)
    # stress_vec_pressure_n1, history_update, internal_energy_n1, temp_k, sdd, junk = evptn_wc.batch_solve(
    #              bdt, drs, svs, vrv, ie, svp, hiv, tk)

    print("Deviatoric Stress + pressure:")
    print(stress_vec_pressure_n1)
    print("Internal energy")
    print(internal_energy_n1)
    print("Temperature (K)")
    print(temp_k)
    print("SDD array")
    print(sdd)
    print("Slip Rates")    
    jax.debug.print("{}", evptn_wc.hist_class.get_slip_rate(history_update))
    print("Deviatoric crystal elastic strain")
    jax.debug.print("{}", evptn_wc.hist_class.get_elas_dev(history_update))
    print("Lattice quaternions")
    jax.debug.print("{}", evptn_wc.hist_class.get_quats(history_update))
    print("Number of function evaluations")
    jax.debug.print("{}", history_update[evptn_wc.hist_class.ind_hist_num_func_evals])
    print("MTan array")
    print(mtan)
