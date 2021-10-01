import ecmech as m
import numpy as np
from scipy.optimize import root

assert m.__version__ == 'dev'

# Helper class for the root solve
class ECMechProb:
    def __init__(self, model_name, var):
        self.myecmech = m.ECMech(model_name, var)

    def getHistInfo(self):
        names, vals, plot, state = self.myecmech.getHistoryInfo()
        return (names, vals, plot, state)
    def setup(self, dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv, sdd, mtanSD):
        self.myecmech.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv, sdd, mtanSD)
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
help(m.ECMech)

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

x = np.zeros(8)

prob.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv, sdd, mtanSD)
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
    # Set this to 0 before each root solve
    x[:] = 0.0
    prob.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv, sdd, mtanSD)
    sol = root(prob.computeRJ, x, jac=True, method='hybr', tol=tolerance)
    eInt, stressSvecP, hist, tkelv, sdd = prob.getState(sol.x, eInt, stressSvecP, hist, tkelv, sdd)