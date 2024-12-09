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

import optimistix as optx

class MemoizeJac:
    """ Decorator that caches the return values of a function returning `(fun, grad)`
        each time it is called. """

    def __init__(self, fun):
        self.fun = fun
        self.jac = None
        self._value = None
        self.x = None

    def _compute_if_needed(self, x, *args):
        if not jnp.all(x == self.x) or self._value is None or self.jac is None:
            self.x = jnp.asarray(x).copy()
            fg = self.fun(x, *args)
            self.jac = fg[1]
            self._value = fg[0]

    def __call__(self, x, *args):
        """ returns the the function value """
        self._compute_if_needed(x, *args)
        return self._value

    def derivative(self, x, *args):
        self._compute_if_needed(x, *args)
        return self.jac

# This is not going to be fun to get into a state where things can be vectorized by JAX...
class SNLSTrDlDenseG:
    class params:
        def __init__(self):
            self.factor = 100.0
            self.maxfev = 200
            self.xtol = jnp.sqrt(jnp.finfo(jnp.float64).eps)

    def __init__(self, functor, xtolerance = jnp.finfo(jnp.float64).eps, ndim = 1, args=()):
        self.fun = MemoizeJac(functor)
        self.jac = self.fun.derivative

        self.nfev  = 0
        self.njev  = 0
        self.iter  = 0
        self.res = 1e20
        self.delta = 1e8
        self.rho_last = 0.0
        self.ndim = ndim
        self.parameters = self.params()
        self.parameters.xtol = xtolerance
        self.success = -10

        self.delta_control = DeltaControl()

        if not isinstance(args, tuple):
            self.args = (args,)
        else:
            self.args = args

    def resetParams(self):
        self.parameters = self.params()
    
    def solveInit(self, x):
        self.delta = self.delta_control.getDeltaInit()
        self.njev = 1
        self.nfev = 1 
        self.success = -10
        self.residual = jnp.copy(self.fun(x, *self.args))
        self.jacobian = jnp.copy(self.jac(x, *self.args))
        self.res = jnp.linalg.norm(self.residual)
        # initialize iteration counter and monitors
        self.iter = 1
        return x

    def solve(self, x):
        x = self.solveInit(x)
        if self.res < self.parameters.xtol:
            self.success = 0
            return (self.success, x)

        res_0 = self.res
        reject_prev = False
        
        nr_step = jnp.zeros(self.ndim)
        grad    = jnp.zeros(self.ndim)
        delta_x = jnp.zeros(self.ndim)
        Jg_2    = 0.0

        for niters in range(self.parameters.maxfev):
            if not reject_prev:
                # This is done outside this step so that these operations can be done with varying solve
                # techniques such as LU/QR or etc...
                grad = self.jacobian.T.dot(self.residual)
                Jg_2 = jnp.dot(self.jacobian.dot(grad), self.jacobian.dot(grad))
                nr_step = jnp.linalg.solve(self.jacobian, self.residual)
                nr_step *= -1.0

            use_nr = False
            # If the step was rejected nrStep will be the same value as previously, and so we can just recalculate nr_norm here.
            nr_norm = jnp.linalg.norm(nr_step)

            # Computes the updated delta x/x, predicated residual error, and whether or not NR method was used.
            use_nr, pred_resid, delta_x, x = self.dogleg(res_0, nr_norm, Jg_2, grad, nr_step, x, use_nr)

            reject_prev = False

            # {
            #    bool rjSuccess = this->computeRJ(residual, Jacobian) ; // at _x
            #    snls::updateDelta<_nDim>(_deltaControl, residual, res_0, pred_resid, nr_norm, _tolerance, use_nr, rjSuccess,
            #                             _delta, _res, _rhoLast, reject_prev, _status, _os);
            #    if(_status != SNLSStatus_t::unConverged) { break; }
            # }
            self.residual = jnp.copy(self.fun(x, *self.args))
            self.jacobian = jnp.copy(self.jac(x, *self.args))
            self.njev += 1
            self.nfev += 1 

            reject_prev = self.update_delta(res_0, pred_resid, nr_norm, use_nr, reject_prev)
            if self.success != -10:
                break

            if reject_prev:
                self.res = res_0
                x -= delta_x

            res_0 = self.res
            
        return (self.success, x)

    def dogleg(self, res_0, nr_norm, Jg_2, grad, nr_step, x, use_nr):
        # No need to do any other calculations if this condition is true
        if ( nr_norm <= self.delta ):
            # use Newton step
            use_nr = True
            delx = nr_step.copy()
            pred_resid = 0.0
        # Find Cauchy point
        else:
            # If we didn't reject things this is the only thing that needs to be updated
            # everything else we should have the info to recompute
            # The nice thing about recomputing is that we can actually define the variables as const
            # to help the compiler out.

            norm2_grad = jnp.dot(grad, grad)
            norm_grad  = jnp.sqrt(norm2_grad)

            alpha = 1.0
            if Jg_2 > 0.0:
                alpha = norm2_grad / Jg_2
            
            norm_grad_inv = 1.0
            if norm_grad > 0.0:
                norm_grad_inv = 1.0 / norm_grad

            norm_s_sd_opt = alpha * norm_grad

            # step along the dogleg path
            if ( norm_s_sd_opt >= self.delta ):
                # use step along steapest descent direction
                delx = -self.delta * norm_grad_inv * grad

                val = -(self.delta * norm_grad) + 0.5 * self.delta * self.delta * Jg_2 * (norm_grad_inv * norm_grad_inv)
                pred_resid = jnp.sqrt(jnp.maximum(2.0 * val + res_0 * res_0, 0.0))
            else:
                qb = 0.0
                qa = 0.0
                p = nr_step + alpha * grad
                qa = jnp.dot(p, p)
                qb = jnp.dot(p, grad)

                # Previously qb = (-p^t g / ||g||) * alpha * ||g|| * 2.0
                # However, we can see that this simplifies a bit and also with the beta term
                # down below we there's a common factor of 2.0 that we can eliminate from everything
                qb *= alpha
                # qc and beta depend on delta
                qc = norm_s_sd_opt * norm_s_sd_opt - self.delta * self.delta
                beta = (qb + jnp.sqrt(qb * qb - qa * qc)) / qa
                beta = jnp.maximum(0.0, jnp.minimum(1.0, beta)) # to deal with any roundoff

                # delx[iX] = alpha*ngrad[iX] + beta*p[iX] = beta*nrStep[iX] - (1.0-beta)*alpha*grad[iX]
                omb  = 1.0 - beta
                omba = omb * alpha
                delx = beta * nr_step - omba * grad
                res_cauchy = res_0
                if Jg_2 > 0.0:
                    res_cauchy = jnp.sqrt(jnp.maximum(0.0, res_0 * res_0 - alpha * norm2_grad))
                pred_resid = omb * res_cauchy

        x += delx
        return (use_nr, pred_resid, delx, x)

    def update_delta(self, res_0, pred_resid, nr_norm, use_nr, reject_prev):
        self.res = jnp.linalg.norm(self.residual)
        # allow to exit now, may have forced one iteration anyway, in which
        # case the delta update can do funny things if the residual was
        # already very small
        if self.res < self.parameters.xtol:
            self.success = 0
            return False
        
        delta_success, reject_prev, self.rho_last, self.delta = self.delta_control.updateDelta(self.delta, self.res, res_0, pred_resid, reject_prev, use_nr, nr_norm, self.rho_last)

        if not delta_success:
            self.success = -20
            return False
        
        return reject_prev

class DeltaControl:
    def __init__(self):
        self.xiLG = 0.75
        self.xiUG = 1.4
        self.xiIncDelta = 1.5
        self.xiLO = 0.35
        self.xiUO = 5.0
        self.xiDecDelta = 0.25
        self.xiForcedIncDelta = 1.2
        self.deltaInit = 1.0
        self.deltaMin = 1e-12
        self.deltaMax = 1e4
        self.rejectResIncrease = True

    def getDeltaInit(self):
        return self.deltaInit

    def decrDelta(self, delta, normfull, took_full):
        success = True

        if took_full:
            delta = jnp.sqrt(delta * self.xiDecDelta * normfull * self.xiDecDelta)
        else:
            delta = delta * self.xiDecDelta

        if delta < self.deltaMin:
            delta = self.deltaMin
            success = False

        return (success, delta)

    def incrDelta(self, delta):
        delta = delta * self.xiIncDelta
        if delta > self.deltaMax:
            delta = deltaMax
        return delta

    def updateDelta(self, delta, res, res_0, pred_res, reject, took_full, normfull, rho):

        success = True
        actual_change = res - res_0
        pred_change = pred_res - res_0
        if pred_change == 0.0:
            if delta >= self.deltaMax:
                # things are going badly enough that the solver should probably stop
                #print("predicted change is zero and delta at max")
                success = False
            else:
                #print("predicted change is zero, forcing delta larger")
                delta = jnp.minimum(delta * self.xiForcedIncDelta, self.deltaMax)
        else:
            rho = actual_change / pred_change

            #print("rho = " + str(rho))

            if (rho > self.xiLG and 
                actual_change < 0.0 and
                rho < self.xiUG
            ):
                if not took_full:
                    #increase delta
                    delta = self.incrDelta(delta)
            elif (rho < self.xiLO or rho > self.xiUO):
                success, delta = self.decrDelta(delta, normfull, took_full)
        reject = False

        if actual_change > 0.0 and self.rejectResIncrease:
            # #print("actual change = " + str(actual_change))
            reject = True

        return (success, reject, rho, delta)

def computeRJ2(x, args=()):
    mlambda = args
    ndim = 8
    r = jnp.zeros(ndim)
    jacob = jnp.zeros((ndim, ndim))
    r = r.at[0].set((3.0 - 2.0 * x[0]) * x[0] - 2.0 * x[1] + 1.0)
    for i in range(1, ndim - 1, 1):
        r = r.at[i].set((3.0 - 2.0 * x[i]) * x[i] - x[i-1] - 2.0 * x[i+1] + 1.0)

    fn = (3.0 - 2.0 * x[-1]) * x[-1] - x[-2] + 1.0
    r = r.at[-1].set((1.0 - mlambda) * fn + mlambda * (fn * fn))

    # F(0) = (3-2*x[0])*x[0] - 2*x[1] + 1
    jacob = jacob.at[0, 0].set(3.0 - 4.0 * x[0])
    jacob = jacob.at[0, 1].set(-2.0)
    # F(i) = (3-2*x[i])*x[i] - x[i-1] - 2*x[i+1] + 1
    for i in range(1, ndim - 1, 1):
        jacob = jacob.at[i, i - 1].set(-1.0)
        jacob = jacob.at[i, i].set(3.0 - 4.0 * x[i])
        jacob = jacob.at[i,i + 1].set(-2.0)

    # F(n-1) = ((3-2*x[n-1])*x[n-1] - x[n-2] + 1)^2;
    fn = (3.0 - 2.0 * x[-1]) * x[-1] - x[-2] + 1.0
    dfndxn = 3.0 - 4.0 * x[-1]
    jacob = jacob.at[-1, -1].set((1.0 - mlambda) * (dfndxn) + mlambda * (2.0 * dfndxn * fn))
    jacob = jacob.at[-1, -2].set((1.0 - mlambda) * (-1.0) + mlambda * (-2.0 * fn))
    
    #print(jnp.linalg.norm(r))
    return (r, jacob)

def computeRJ3(x, args):
    mlambda = args
    ndim = 8
    r = jnp.zeros(ndim)
    r = r.at[0].set((3.0 - 2.0 * x[0]) * x[0] - 2.0 * x[1] + 1.0)
    for i in range(1, ndim - 1, 1):
        r = r.at[i].set((3.0 - 2.0 * x[i]) * x[i] - x[i-1] - 2.0 * x[i+1] + 1.0)

    fn1 = (3.0 - 2.0 * x[-1]) * x[-1] - x[-2] + 1.0
    r = r.at[-1].set((1.0 - mlambda) * fn1 + mlambda * (fn1 * fn1))

    return r

if __name__ == "__main__":
    x = jnp.ones(8) * 0.0
    args = (0.99999999)

    solver = optx.Dogleg(rtol=1e-6, atol=1e-8)
    sol = optx.root_find(computeRJ3, solver=solver, y0=x, args=args)
    print(sol.stats)
    xs = sol.value

    @partial(jax.jit, static_argnums=(1,2))
    def root_find(x, fn, args):
        solver = optx.Dogleg(rtol=1e-6, atol=1e-8)
        return optx.root_find(fn=fn, solver=solver, y0=x, args=args).value

    xs2 = root_find(x, computeRJ3, args)

    print(xs)
    print(xs2)
    print(np.linalg.norm(computeRJ3(xs, args)))

    x = jnp.ones(8) * 0.0
    solver1 = SNLSTrDlDenseG(computeRJ2, xtolerance=1e-12, ndim=x.shape[0], args=args)
    solver1.delta_control.deltaInit = 100.0
    status, xs = solver1.solve(x)
    print(status, xs)
    res, _ = computeRJ2(xs, args)
    print(np.linalg.norm(res))
    print(solver1.nfev, solver1.njev)
    print()
