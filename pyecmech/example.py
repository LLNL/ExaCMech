"""
Minimal end-to-end usage example for pyecmech.

Builds a linear-Voce FCC model ("voce_fcc_norm") with example OFHC copper parameters,
drives it through a simple monotonic uniaxial-tension deformation history over 41
explicit time steps, and prints the resulting axial Cauchy stress at each step.

See ECMechProb.solve() below for the array shapes pyECMech.solve() expects, and
ecmechpy.hpp / ecmech_pybind11.cpp for the full documentation of the underlying C++
binding.
"""

import pyecmech as m
import numpy as np

assert m.__version__ == 'dev'


# Thin convenience wrapper pairing a pyecmech.pyECMech instance with the dimension
# constants needed to correctly shape the numpy arrays passed to solve().
class ECMechProb:
    def __init__(self, model_name, var):
        self.myecmech = m.pyECMech(model_name, var)
        # Find the number of history variables and store that for future use
        self.nhist = self.myecmech.getNumberHistory()
        self.ntvec = m.constants.ntvec
        self.nsvec = m.constants.nsvec
        self.nsvec2 = m.constants.nsvec2
        self.nvr = m.constants.nvr
        self.ne = m.constants.ne
        self.nsvp = m.constants.nsvp
        self.nwvec = m.constants.nwvec
        self.nsdd = m.constants.nsdd

    def getHistInfo(self):
        names, vals, plot, state = self.myecmech.getHistoryInfo()
        return (names, vals, plot, state)

    def solve(self, dt, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k):
        '''
            Solve does a per time step solve of the material update for all points inputted.
            A few things to note:
            def_rate_dev6_vol_sample has dimensions npts x self.nsvp input
            spin_vec_sample has dimesnions npts x self.nwvec input
            volRatio has dimensions npts x self.nvr input
            internal_energy has dimensions npts x self.ne input/output
            cauchy_stress_dev6_pressure has dimensions npts x self.nsvec input/output
            hist has dimensions npts x self.nhist input/output
            temp_k has dimensions npts x 1 input/output
            sdd has dimensions npts x self.nsdd output

            If you pass in 1D arrays we will promote them to 2D arrays.
        '''

        def_rate_dev6_vol_sample = np.atleast_2d(def_rate_dev6_vol_sample)
        spin_vec_sample = np.atleast_2d(spin_vec_sample)
        volRatio = np.atleast_2d(volRatio)
        internal_energy = np.atleast_2d(internal_energy)
        cauchy_stress_dev6_pressure = np.atleast_2d(cauchy_stress_dev6_pressure)
        hist = np.atleast_2d(hist)
        temp_k = np.atleast_2d(temp_k)

        npts = internal_energy.shape[0]
        sdd = np.zeros((npts, self.nsdd))
        self.myecmech.solve(dt, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd, npts)
        return (internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd)

# Prints out function documentation of the module
# help(m)

# OFHC copper parameters for "voce_fcc_norm" ("evptn_FCC_A": FCC slip geometry (0 of its
# own params) + linear-Voce power-law kinetics + cubic elastic constants + simple EOS).
# Taken from parameters in the ExaConstit test suite. Order below follows
# matModel::initFromParams() in ECMech_evptnWrap.h: density0, cvav, tolerance, then
# SlipGeom params (none for FCC), then ThermoElastN, Kinetics, and finally the remaining
# EOS params -- see ECMech_eosSimple.h / ECMech_elastic.h / kinetics/ECMech_kinetics_VocePL.h
# for what each model's own slice means.
var = np.asarray([
8.920e-6,    # rho0    -- reference density
0.003435984, # cvav    -- specific heat
1.0e-10,     # tolerance -- solver tolerance
168.4e0,     # C11     -- cubic elastic constant
121.4e0,     # C12     -- cubic elastic constant
75.2e0,      # C44     -- cubic elastic constant
44.0e0,      # mu      -- shear modulus used by the power-law slip kinetics
0.02e0,      # xm      -- rate-sensitivity exponent (power-law slip kinetics)
1.0e0,       # gam_w   -- reference/normalizing shear rate (power-law slip kinetics)
400.0e-3,    # h0      -- initial (linear) Voce hardening rate
17.0e-3,     # tausi   -- initial CRSS (Voce hardening)
122.4e-3,    # taus0   -- Voce saturation-stress reference value
0.0,         # xms     -- Voce saturation-stress rate-sensitivity exponent
5.0e9,       # gamss0  -- reference shear rate for the saturation stress
17.0e-3,     # hdn_init -- initial hardening state (matches tausi here)
0.0,         # Gamma   -- EOS Gruneisen parameter (remaining EOS params; rho0/K0/cvav are
             #            supplied to the EOS internally from the values above)
-1.0307952   # e0      -- EOS reference/offset energy
])

prob = ECMechProb("voce_fcc_norm", var)

# # Our various input parameters
dt = 0.1
tolerance = 1e-10
def_rate_dev6_vol_sample = np.zeros(7)
# Just a simple monotonic tension example in the x direction. def_rate_dev6_vol_sample is
# laid out as [deviatoric 6-vector, volumetric rate] (pyecmech.constants.nsvp == 7 wide;
# see matModelBase::getResponseECM's def_rate_d6vV doc). Indices 0-2 encode the traceless
# (deviatoric) part of a uniaxial stretch along x (-1/3 on the two transverse directions
# for every +1 on x), while index 6 carries the accompanying volumetric strain-rate
# component (nonzero here, so volRatio below does drift away from 1.0 over the run); both
# are then scaled down to a small strain rate.
d_tr = 1.0 / 3.0
def_rate_dev6_vol_sample[0] = 1.0 - d_tr
def_rate_dev6_vol_sample[1] = -d_tr
def_rate_dev6_vol_sample[2] = -d_tr
def_rate_dev6_vol_sample[6] = 3.0 * d_tr
def_rate_dev6_vol_sample[:] *= 0.001

cauchy_stress_dev6_pressure = np.zeros(7)
# This would control the spin of the problem if we wanted to
spin_vec_sample = np.zeros(3)
internal_energy = np.zeros(1)
# [rel_vol_n, rel_vol_n+1, rate, delta] -- see ECMech_const.h's `nvr` doc. Both start at
# 1.0 (undeformed reference state); the loop below advances rel_vol_n+1 each step by
# integrating the volumetric-rate component (index 6) of def_rate_dev6_vol_sample.
volRatio = np.asarray([1.0, 1.0, 0.0, 0.0])

temp_k = 300.
sdd = np.asarray([0, 0])
mtanSD = np.zeros(36)

histNames, histVals, histPlot, histState = prob.getHistInfo()

# Just so we can see what the history names and values are
# print(histNames)
# print(histVals)

hist = np.copy(histVals)
# How to iterate over multiple time steps
for i in range(41):
    # This is pulled from how the test_px does things. Roll last step's rel_vol_n+1 into
    # this step's rel_vol_n, integrate the (constant) volumetric strain rate to get the
    # new rel_vol_n+1, then recompute the derived rate/delta slots -- see volRatio's
    # [rel_vol_n, rel_vol_n+1, rate, delta] layout noted above.
    volRatio[0] = volRatio[1]
    volRatio[1] = volRatio[0] * np.exp(def_rate_dev6_vol_sample[6] * dt)
    volRatio[3] = volRatio[1] - volRatio[0]
    volRatio[2] = volRatio[3] / (dt * 0.5 * (volRatio[0] + volRatio[1]))

    # An example of using prob.solve()
    internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd = prob.solve(dt, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k)

    print(cauchy_stress_dev6_pressure[:,0])