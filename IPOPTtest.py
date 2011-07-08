'''
Solves problem number 71 from the Hock-Schittkowsky test suite.

    min     x1*x4*(x1 + x2 + x3) + x3
    s.t.:   
            x1**2 + x2**2 + x3**2 + x4**2 - 40 = 0
            x1*x2*x3*x4 - 25 >= 0
            0 <= xi <= 5,  i = 1,2,3,4
    

    x0 = [1, 5, 5, 1]    (Feasible)

    f* = 17.0140173
    x* = [1, 4.7429994, 3.8211503, 1.3794082]
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time
import pdb
import numpy
# =============================================================================
# Extension modules
# =============================================================================
#from pyOpt import *
from pyOpt import Optimization
from pyOpt import IPOPT


# =============================================================================
# 
# =============================================================================
def objfunc(x):
    
    f = x[0]*x[3]*(x[0] + x[1] + x[2]) + x[2]
    g = [0.0]*2
    g[0] = x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2
    g[1] = x[0]*x[1]*x[2]*x[3]
    
    fail = 0
    return f,g, fail
    

def gradfunc(x,f,g):
    
    g_obj = [0.0]*4
    g_obj[0] = x[3]*(2.*x[0] + x[1] + x[2])
    g_obj[1] = x[0]*x[3]
    g_obj[2] = x[0]*x[3] + 1
    g_obj[3] = x[0]*(x[0] + x[1] + x[2])
    
    g_con = numpy.zeros([2,4])
    g_con[0][0] = 2.*x[0]
    g_con[0][1] = 2.*x[1]
    g_con[0][2] = 2.*x[2]
    g_con[0][3] = 2.*x[3]
    g_con[1][0] = x[1]*x[2]*x[3]
    g_con[1][1] = x[0]*x[2]*x[3]
    g_con[1][2] = x[0]*x[1]*x[3]
    g_con[1][3] = x[0]*x[1]*x[2]
   
    fail = 0
    
    return g_obj,g_con,fail
    

# =============================================================================
# 
# =============================================================================

# Instantiate Optimization Problem
opt_prob = Optimization('HS 071',objfunc)
opt_prob.addVar('x1','c',value=1.,lower=1.,upper=5.)
opt_prob.addVar('x2','c',value=5.,lower=1.,upper=5.)
opt_prob.addVar('x3','c',value=5.,lower=1.,upper=5.)
opt_prob.addVar('x4','c',value=1.,lower=1.,upper=5.)
opt_prob.addObj('f')
opt_prob.addCon('g1','e',equal=40.)
opt_prob.addCon('g2','i',lower=25.,upper=numpy.inf)

print opt_prob

# Instantiate Optimizer (IPOPT)
ipopt = IPOPT()

# Solve Problem with Optimizer Using Finite Differences
ipopt.setOption('IFILE','ipopt.out')
ipopt.setOption('max_iter',2)
ipopt.setOption('linear_system_scaling','none')
ipopt.setOption('tol',1e-2)

ipopt(opt_prob,sens_type='FD',store_hst=True)
print opt_prob.solution(0)

