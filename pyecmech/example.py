import pyecmech as m
import numpy as np

assert m.__version__ == 'dev'


# Helper class for things
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

prob = ECMechProb("voce_fcc_norm", var)

# # Our various input parameters
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
    # This is pulled from how the test_px does things
    volRatio[0] = volRatio[1]
    volRatio[1] = volRatio[0] * np.exp(def_rate_dev6_vol_sample[6] * dt)
    volRatio[3] = volRatio[1] - volRatio[0]
    volRatio[2] = volRatio[3] / (dt * 0.5 * (volRatio[0] + volRatio[1]))

    # An example of using prob.solve()
    internal_energy, cauchy_stress_dev6_pressure, hist, temp_k, sdd = prob.solve(dt, def_rate_dev6_vol_sample, spin_vec_sample, volRatio, internal_energy, cauchy_stress_dev6_pressure, hist, temp_k)

    print(cauchy_stress_dev6_pressure[:,0])