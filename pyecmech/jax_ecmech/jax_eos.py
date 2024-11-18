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

class eosSimple:
    def __init__(
                 self,
                 parameters,
                 is_isothermal
                 ):
        self.num_params = 5
        self.num_hist = 0

        self.isothermal = is_isothermal

        self.density_0 = parameters["init_density"]
        self.bulk_modulus = parameters["bulk_modulus_0"]
        self.cvav = parameters["cvav"]
        self.gamma = parameters["eos_gamma"]
        self.cold_energy_0 = parameters["eos_cold_energy_0"]

        self.dtde = 1.0 / self.cvav
        self.temp_k_init = -self.cold_energy_0 * self.dtde

    def get_parameters(self, parameters):
        parameters["init_density"] = self.density_init
        parameters["bulk_modulus_0"] = self.bulk_modulus
        parameters["cvav"] = self.cvav
        parameters["eos_gamma"] = self.gamma
        parameters["eos_cold_energy_0"] = self.cold_energy_0

        return parameters

    def eval_pressure_temp(self, volume, energy):
        mu = 1.0 / volume - 1.0

        pressure = self.bulk_modulus * mu
        temp_k = self.temp_k_init

        pressure, temp_k = jax.lax.cond(
            not self.isothermal,
            lambda: (pressure + self.gamma * energy, temp_k + self.dtde * energy),
            lambda: (pressure, temp_k)
        )

        # if not self.isothermal:
        #     pressure += self.gamma * energy
        #     temp_k += self.dtde * energy
        return (pressure, temp_k)

    def eval_temp(self, energy):
        temp_k = self.temp_k_init

        temp_k = jax.lax.cond(
            not self.isothermal,
            lambda: temp_k + self.dtde * energy,
            lambda: temp_k
        )

        # if not self.isothermal:
        #     temp_k += self.dtde * energy
        return temp_k

    def eval_pressure_temp_diff(self, volume, energy):
        eta = 1.0 / volume
        mu = eta - 1.0

        temp_k = self.eval_temp(energy)
        bulk_mod_new = self.bulk_modulus * eta

        pressure = self.bulk_modulus * mu
        dtde = self.dtde
        dpde = 0.0

        dtde, dpde, pressure = jax.lax.cond(
            self.isothermal,
            lambda: (dtde * 1e-8, dpde, pressure),
            lambda: (dtde, self.gamma, pressure + self.gamma * energy)
        )
        
        # if self.isothermal:
        #     dtde *= 1e-8
        # else:
        #     pressure += self.gamma * energy
        #     dpde = self.gamma

        return (pressure, temp_k, bulk_mod_new, dpde, dtde)

def update_simple(eos_class, vol_new, vol_incr, energy_old, pressure_old):
    energy_new = energy_old - vol_incr * pressure_old
    pressure, temp_k, bulk_mod_new, dpde, dtde = eos_class.eval_pressure_temp_diff(vol_new, energy_new)

    dpdv = -bulk_mod_new / vol_new

    bulk_mod_new = jnp.max(jnp.asarray([1e-5 * eos_class.bulk_modulus, bulk_mod_new + dpde * pressure_old * vol_new]))

    return (temp_k, pressure, energy_new, bulk_mod_new)