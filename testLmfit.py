#!/usr/bin/env python
#<examples/doc_withreport.py>

#==============================================================================
# from __future__ import print_function
# from lmfit import Parameters, minimize, fit_report
# from numpy import random, linspace, pi, exp, sin, sign
# 
# 
# p_true = Parameters()
# p_true.add('amp', value=14.0)
# p_true.add('period', value=5.46)
# p_true.add('shift', value=0.123)
# p_true.add('decay', value=0.032)
# 
# def residual(pars, x, data=None):
#     vals = pars.valuesdict()
#     amp =  vals['amp']
#     per =  vals['period']
#     shift = vals['shift']
#     decay = vals['decay']
# 
#     if abs(shift) > pi/2:
#         shift = shift - sign(shift)*pi
#     model = amp * sin(shift + x/per) * exp(-x*x*decay*decay)
#     if data is None:
#         return model
#     return (model - data)
# 
# n = 1001
# xmin = 0.
# xmax = 250.0
# 
# random.seed(0)
# 
# noise = random.normal(scale=0.7215, size=n)
# x     = linspace(xmin, xmax, n)
# data  = residual(p_true, x) + noise
# 
# fit_params = Parameters()
# fit_params.add('amp', value=13.0)
# fit_params.add('period', value=2)
# fit_params.add('shift', value=0.0)
# fit_params.add('decay', value=0.02)
# 
# out = minimize(residual, fit_params, args=(x,), kws={'data':data})
# print(fit_report(out))


# #<end of examples/doc_withreport.py>

#==============================================================================



#<examples/doc_basic.py>
from lmfit import minimize, Parameters, Parameter, report_fit
import numpy as np

# create data to be fitted
x = np.linspace(0, 15, 301)
data = (5. * np.sin(2 * x - 0.1) * np.exp(-x*x*0.025) +
        np.random.normal(size=len(x), scale=0.2) )

# define objective function: returns the array to be minimized
def fcn2min(params, x, data):
    """ model decaying sine wave, subtract data"""
    amp = params['amp'].value
    shift = params['shift'].value
    omega = params['omega'].value
    decay = params['decay'].value

    model = amp * np.sin(x * omega + shift) * np.exp(-x*x*decay)
    return model - data

# create a set of Parameters
params = Parameters()
params.add('amp',   value= 10,  min=0)
params.add('decay', value= 0.1)
params.add('shift', value= 0.0, min=-np.pi/2., max=np.pi/2)
params.add('omega', value= 3.0)


# do fit, here with leastsq model
result = minimize(fcn2min, params, args=(x, data))

# calculate final result
final = data + result.residual

# write error report
report_fit(result.params)

# try to plot results
try:
    import pylab
    pylab.plot(x, data, 'k+')
    pylab.plot(x, final, 'r')
    pylab.show()
except:
    pass

#<end of examples/doc_basic.py>