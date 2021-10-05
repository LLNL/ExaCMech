import ecmech as m
import numpy as np
from scipy.optimize import root

assert m.__version__ == 'dev'

# Helper class for the root solve
class ECMechProb:
    def __init__(self, model_name, var):
        self.myecmech = m.ECMech(model_name, var)
        # Find the number of history variables and store that for future use
        names, vals, plot, state = self.myecmech.getHistoryInfo()
        self.nhist = len(names)
        # Maybe have these provided by myecmech?
        self.ntvec = 5
        self.nsvec = 6
        self.nsvec2 = 36
        self.nvr = 4
        self.ne = 1
        self.nsvp = 7
        self.nwvec = 3
        self.nsdd = 2

    def solve(self, dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv):
        '''
            Solve does a per time step solve of the material update for all points inputted.
            A few things to note:
            d_svec_kk_sm has dimensions npts x self.nsvp input
            w_veccp_sm has dimesnions npts x self.nwvec input
            volRatio has dimensions npts x self.nvr input
            eInt has dimensions npts x self.ne input/output
            stressSvecP has dimensions npts x self.nsvec input/output
            hist has dimensions npts x self.nhist input/output
            tkelv has dimensions npts x 1 input/output
            sdd has dimensions npts x self.nsdd output

            If you pass in 1D arrays we will promote them to 2D arrays.
        '''
        d_svec_kk_sm = np.atleast_2d(d_svec_kk_sm)
        w_veccp_sm = np.atleast_2d(w_veccp_sm)
        volRatio = np.atleast_2d(volRatio)
        eInt = np.atleast_2d(eInt)
        stressSvecP = np.atleast_2d(stressSvecP)
        hist = np.atleast_2d(hist)
        tkelv = np.atleast_2d(tkelv)

        npts = eInt.shape[0]
        sdd = np.zeros((npts, self.nsdd))

        x = np.zeros(8)

        for i in range(npts):
            x[:] = 0.0
            self.setup(dt, tolerance, np.squeeze(d_svec_kk_sm[i, :]), np.squeeze(w_veccp_sm[i, :]), np.squeeze(volRatio[i, :]), np.squeeze(eInt[i, :]), np.squeeze(stressSvecP[i, :]), np.squeeze(hist[i, :]), np.squeeze(tkelv[i, :]))
            sol = root(self.computeRJ, x, jac=True, method='hybr', tol=tolerance)
            # If you want to check the success of the solver you can find that using
            # sol.success
            eInt[i, :], stressSvecP[i, :], hist[i, :], tkelv[i, :], sdd[i, :] = prob.getState(sol.x,  np.squeeze(eInt[i, :]), np.squeeze(stressSvecP[i, :]), np.squeeze(hist[i, :]), np.squeeze(tkelv[i, :]), np.squeeze(sdd[i, :]))

        return (eInt, stressSvecP, hist, tkelv, sdd)

    def getHistInfo(self):
        names, vals, plot, state = self.myecmech.getHistoryInfo()
        return (names, vals, plot, state)
    def setup(self, dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv):
        self.myecmech.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv)
    def computeRJ(self, x):

        ndim = x.shape[0]
        resid = np.zeros(ndim)
        jacob = np.zeros((ndim,ndim))
        self.myecmech.computeRJ(resid, jacob, x)

        return (resid, jacob)

    def getState(self, x, eInt, stressSvecP, hist, tkelv, sdd):
        self.myecmech.getState(x, eInt, stressSvecP, hist, tkelv, sdd)
        return (eInt, stressSvecP, hist, tkelv, sdd)

# Prints out function documentation of the module
# help(m.ECMech)

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

# Our various input parameters
dt = 0.1
tolerance = 1e-10
d_svec_kk_sm = np.zeros(7)
# Just a simple monotonic tension example in the x direction
d_tr = 1.0 / 3.0
d_svec_kk_sm[0] = 1.0 - d_tr
d_svec_kk_sm[1] = -d_tr
d_svec_kk_sm[2] = -d_tr
d_svec_kk_sm[6] = 3.0 * d_tr
d_svec_kk_sm[:] *= 0.001

stressSvecP = np.zeros(7)
# This would control the spin of the problem if we wanted to 
w_veccp_sm = np.zeros(3)
eInt = np.zeros(1)
volRatio = np.asarray([1.0, 1.0, 0.0, 0.0])

volRatio[0] = volRatio[1]
volRatio[1] = volRatio[0] * np.exp(d_svec_kk_sm[6] * dt)
volRatio[3] = volRatio[1] - volRatio[0]
volRatio[2] = volRatio[3] / (dt * 0.5 * (volRatio[0] + volRatio[1]))

tkelv = 300.
sdd = np.asarray([0, 0])
mtanSD = np.zeros(36)

histNames, histVals, histPlot, histState = prob.getHistInfo()

# Just so we can see what the history names are
print(histNames)

hist = np.copy(histVals)
hist_old = np.copy(hist)

# An example of how to manually solve for things if you want to play around with different
# solver options or if you just don't want to use the ECMechProb.solve() function
x = np.zeros(8)
prob.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv)
sol = root(prob.computeRJ, x, jac=True, method='hybr', tol=tolerance)
# If you want to check the success of the solver you can find that using
# sol.success
eInt, stressSvecP, hist, tkelv, sdd = prob.getState(sol.x, eInt, stressSvecP, hist, tkelv, sdd)

# print(hist)
# How to iterate over multiple time steps
for i in range(40):
    # This is pulled from how the test_px does things
    volRatio[0] = volRatio[1]
    volRatio[1] = volRatio[0] * np.exp(d_svec_kk_sm[6] * dt)
    volRatio[3] = volRatio[1] - volRatio[0]
    volRatio[2] = volRatio[3] / (dt * 0.5 * (volRatio[0] + volRatio[1]))
    # An example of using the prob.solve() version of things rather than
    # doing it by hand
    eInt, stressSvecP, hist, tkelv, sdd = prob.solve(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv)

    print(stressSvecP)