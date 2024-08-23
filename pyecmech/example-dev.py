import pyecmech as m
import numpy as np
from scipy.optimize import root

assert m.__version__ == 'dev'

# This file shows more how to use a developer set of tools to explore
# the effects of different solvers or other features related to the
# solve of a material point
# It will be enabled by a cmake flag but currently that has not been added

# Helper class for the root solve
class ECMechProbDev:
    def __init__(self, model_name, var):
        self.myecmech = m.pyECMechDev(model_name, var)
        # Find the number of history variables and store that for future use
        names, vals, plot, state = self.myecmech.getHistoryInfo()
        self.nhist = len(names)
        self.ntvec = m.constants.ntvec
        self.nsvec = m.constants.nsvec
        self.nsvec2 = m.constants.nsvec2
        self.nvr = m.constants.nvr
        self.ne = m.constants.ne
        self.nsvp = m.constants.nsvp
        self.nwvec = m.constants.nwvec
        self.nsdd = m.constants.nsdd

    def solve(self, dt, tolerance, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k):
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

        x = np.zeros(8)

        for i in range(npts):
            x[:] = 0.0
            self.setup(dt, tolerance, np.squeeze(def_rate_dev6_vol_sample[i, :]), np.squeeze(spin_vec_sample[i, :]), np.squeeze(volRatio[i, :]), np.squeeze(internal_energy[i, :]), np.squeeze(cauchy_stress_dev6_pressure[i, :]), np.squeeze(hist[i, :]), np.squeeze(temp_k[i, :]))
            sol = root(self.computeRJ, x, jac=True, method='hybr', tol=tolerance)
            # If you want to check the success of the solver you can find that using
            # sol.success
            internal_energy[i, :], cauchy_stress_dev6_pressure[i, :], hist[i, :], temp_k[i, :], sdd[i, :] = prob.getState(sol.x,  np.squeeze(internal_energy[i, :]), np.squeeze(cauchy_stress_dev6_pressure[i, :]), np.squeeze(hist[i, :]), np.squeeze(temp_k[i, :]), np.squeeze(sdd[i, :]))

        return (internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd)

    def getHistInfo(self):
        names, vals, plot, state = self.myecmech.getHistoryInfo()
        return (names, vals, plot, state)
    def setup(self, dt, tolerance, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k):
        self.myecmech.setup(dt, tolerance, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k)
    def computeRJ(self, x):

        ndim = x.shape[0]
        resid = np.zeros(ndim)
        jacob = np.zeros((ndim,ndim))
        self.myecmech.computeRJ(resid, jacob, x)

        return (resid, jacob)

    def getState(self, x, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd):
        self.myecmech.getState(x, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd)
        return (internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd)

# Prints out function documentation of the module
# help(m)

# OFHC copper parameters
# Taken from parameters in the ExaConstit test suite
var = np.asarray([
8.920e-6,
0.003435984,
1.0e-10,
168.4e0,
121.4e0,
75.2e0,
44.0e0,
0.02e0,
1.0e0,
400.0e-3,
17.0e-3,
122.4e-3,
0.0,
5.0e9,
17.0e-3,
0.0,
-1.0307952
])

prob = ECMechProbDev("voce_fcc_norm", var)

# Our various input parameters
dt = 0.1
tolerance = 1e-10
def_rate_dev6_vol_sample = np.zeros(7)
# Just a simple monotonic tension example in the x direction
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
volRatio = np.asarray([1.0, 1.0, 0.0, 0.0])

volRatio[0] = volRatio[1]
volRatio[1] = volRatio[0] * np.exp(def_rate_dev6_vol_sample[6] * dt)
volRatio[3] = volRatio[1] - volRatio[0]
volRatio[2] = volRatio[3] / (dt * 0.5 * (volRatio[0] + volRatio[1]))

temp_k = 300.
sdd = np.asarray([0, 0])
mtanSD = np.zeros(36)

histNames, histVals, histPlot, histState = prob.getHistInfo()

# Just so we can see what the history names are
# print(histNames)

hist = np.copy(histVals)
hist_old = np.copy(hist)

# An example of how to manually solve for things if you want to play around with different
# solver options or if you just don't want to use the ECMechProb.solve() function
x = np.zeros(8)
prob.setup(dt, tolerance, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k)
sol = root(prob.computeRJ, x, jac=True, method='hybr', tol=tolerance)
# If you want to check the success of the solver you can find that using
# sol.success
internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd = prob.getState(sol.x, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd)
print(cauchy_stress_dev6_pressure[0])
# print(hist)
# How to iterate over multiple time steps
for i in range(40):
    # This is pulled from how the test_px does things
    volRatio[0] = volRatio[1]
    volRatio[1] = volRatio[0] * np.exp(def_rate_dev6_vol_sample[6] * dt)
    volRatio[3] = volRatio[1] - volRatio[0]
    volRatio[2] = volRatio[3] / (dt * 0.5 * (volRatio[0] + volRatio[1]))
    # An example of using the prob.solve() version of things rather than
    # doing it by hand
    internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd = prob.solve(dt, tolerance, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k)

    print(cauchy_stress_dev6_pressure[:,0])