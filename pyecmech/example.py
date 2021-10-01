import ecmech as m
import numpy as np
from scipy.optimize import root

assert m.__version__ == 'dev'

var = np.asarray([])

dt = 0.2
d_svec_kk_sm = np.asarray([-0.0903958, 0.0275513, 0.0628445, 0.0277765, 0.0120316, 0.00948184, -0.00241654])
stressSvecP = np.zeros(7) 
w_veccp_sm = np.zeros(3)
eInt = np.zeros(1)
volRatio = np.asarray([1.0, 1.0, 0.0, 0.0])
hist = np.asarray([])
hist_old = np.copy(hist)

x = np.zeros(8)

prob = ECMechProb("voce_fcc_norm", var)
prob.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv, sdd, mtanSD)

sol = root(prob.computeRJ, x, jac=True, method='hybr', tol=1e-10)
eInt, stressSvecP, hist, tkelv, sdd = prob.getState(sol.x, eInt, stressSvecP, hist, tkelv, sdd)

# print(hist)
# How to iterate over multiple time steps
for i in range(20):
    prob.setup(dt, tolerance, d_svec_kk_sm, w_veccp_sm, volRatio, eInt, stressSvecP, hist, tkelv, sdd, mtanSD)
    sol = root(prob.computeRJ, x, jac=True, method='hybr', tol=1e-10)
    eInt, stressSvecP, hist, tkelv, sdd = prob.getState(sol.x, eInt, stressSvecP, hist, tkelv, sdd)