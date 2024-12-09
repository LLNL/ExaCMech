#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 09:38:58 2023

@author: carson16
"""

import numpy as np

import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

try:
    import pyecmech as pecm
    NWVEC = pecm.constants.nwvec
    NTVEC = pecm.constants.ntvec
    NSVP = pecm.constants.nsvp
    NSDD = pecm.constants.nsdd
    QDIM = pecm.constants.qdim
    DBL_TINY_SQRT = pecm.constants.dbl_tiny_sqrt
    GAM_RATIO_OVFFX = pecm.constants.gam_ratio_ovffx
    GAM_RATIO_MIN = pecm.constants.gam_ratio_min
    GAM_RATIO_OVF = pecm.constants.gam_ratio_ovf
    LN_GAM_RATIO_MIN = pecm.constants.ln_gam_ratio_min
    ELAS_SCALE = pecm.constants.e_scale
    ROT_SCALE = pecm.constants.r_scale
except:
    NWVEC = int(3)
    NTVEC = int(5)
    NSVP = int(7)
    NSDD = int(2)
    QDIM = int(4)
    DBL_TINY_SQRT = 1.0e-90
    GAM_RATIO_OVFFX = 1.0e45
    GAM_RATIO_MIN = 1.0e-60
    GAM_RATIO_OVF = 1.0e60
    LN_GAM_RATIO_MIN = -138.15
    ELAS_SCALE = 5e-4
    ROT_SCALE = 0.01
GAM_RATIO_MAX = 1.0e30
EPS = 2.22e-16

GAM_RATIO_MAX = 1.0e30
EPS = 2.22e-16

class HistClass:
    def __init__(self,
                 slip_geom_class, 
                 slip_kinetics_class,
                 thermo_elas_class,
                 eos_class):

        self.num_hist_auxillary = 4

        self.ind_hist_lba = 0
        self.ind_hist_shear_rate_eff = self.ind_hist_lba + 0
        self.ind_hist_shear_eff = self.ind_hist_shear_rate_eff + 1
        self.ind_flow_strength = self.ind_hist_lba + 2
        self.ind_hist_num_func_evals = self.ind_hist_lba + 3

        self.ind_hist_elas_strain = self.num_hist_auxillary
        self.ind_hist_elas_strain_end = self.ind_hist_elas_strain + pecm.constants.ntvec

        self.ind_hist_quats = self.num_hist_auxillary + pecm.constants.ntvec
        self.ind_hist_quats_end = self.ind_hist_quats + pecm.constants.qdim

        self.ind_hist_hard = self.num_hist_auxillary + pecm.constants.ntvec + pecm.constants.qdim
        self.ind_hist_hard_end = self.ind_hist_hard + slip_kinetics_class.num_hard

        self.ind_hist_slip = self.ind_hist_hard + slip_kinetics_class.num_hard
        self.ind_hist_slip_end = self.ind_hist_slip + slip_geom_class.num_slip_systems

        self.num_hist = self.ind_hist_slip + slip_kinetics_class.num_hard + slip_geom_class.num_slip_systems
        

    def get_hard_state(self, hist):
        ibeg = self.ind_hist_hard
        iend = self.ind_hist_hard_end
        return hist[ibeg:iend]

    def get_elas_dev(self, hist):
        ibeg = self.ind_hist_elas_strain
        iend = self.ind_hist_elas_strain_end
        return hist[ibeg:iend]

    def get_slip_rate(self, hist):
        ibeg = self.ind_hist_slip
        iend = self.ind_hist_slip_end
        return hist[ibeg:iend]    

    def get_quats(self, hist):
        ibeg = self.ind_hist_quats
        iend = self.ind_hist_quats_end
        return hist[ibeg:iend]

    def get_shear_eff(self, hist):
        return hist[self.ind_hist_shear_eff]

    def pack_history_vars(self, elas_dev, quats, hard_state, slip_rate, shear_rate_eff, shear_eff, flow_strength, solver_iters):
        return jnp.hstack((shear_rate_eff, shear_eff, flow_strength, solver_iters, elas_dev, quats, hard_state, slip_rate))

