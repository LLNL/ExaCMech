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


class thermoElastCubic:
    def __init__(
                 self,
                 parameters
                 ):
        
        self.num_params = 3
        self.c11 = parameters["C11"]
        self.c12 = parameters["C12"]
        self.c44 = parameters["C44"]

        self.c_diags = jnp.asarray([
                            self.c11 - self.c12,
                            self.c11 - self.c12,
                            2.0 * self.c44,
                            2.0 * self.c44,
                            2.0 * self.c44
                        ])

        self.bulk_mod = 1.0 / 3.0 * (self.c11 + 2.0 * self.c12)
        # average of c_diag terms
        self.shear_mod = (2.0 * self.c11 - 2.0 * self.c12 + 6.0 * self.c44) * 0.2

    def get_parameters(self, parameters):
        parameters["C11"] = self.c11
        parameters["C12"] = self.c12
        parameters["C44"] = self.c44

        return parameters

    def eval(self, elas_dev_press_vec, pressure, eVref):
        jacob = jnp.exp(jnp.sqrt(3.0) * elas_dev_press_vec[-1])
        kirchoff_bulk = -jnp.sqrt(3.0) * jacob * pressure
        kirchoff_dev = self.c_diags * elas_dev_press_vec[0:-1]
        return jnp.hstack((kirchoff_dev, kirchoff_bulk))

    def calc_dtau_depsilon(dtau_deps_mat, schmid_system_p_vecs, inv_a_vol):
        return dtau_deps_mat + (schmid_system_p_vecs * self.c_diags) * inv_a_vol

    def get_cauchy_stress(kirchoff_stress, inv_det_vol):
        return inv_det_vol * kirchoff_stress
    
    def mult_cauchy_deriv(gen_dev_matrix, inv_det_vol, inv_a_vol):
        dev_mat_contrib = (gen_dev_matrix * self.c_diags) * inv_det_vol * inv_a_vol
        #Pad the pressure related terms with zeros as cubic materials don't contribute
        #anything to those areas
        return jnp.vstack((jnp.c_[dev_mat_contrib, jnp.zeros(5)], jnp.zeros(6)))