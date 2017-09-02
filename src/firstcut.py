"""First cut"""
from numpy import *
# @timer
def firstcut():
    B, r, w = 0.95, 0.04, 1
    p = [0.1, 0.9]

    m = 10
    k = linspace(0, 5, m)
    kp = k 

    J_old = zeros((2, m))
    J_new = zeros((2, m))
    pol = array([zeros((3, m)), zeros((3, m))])

    while True:
        J_old = J_new.copy()
        for i in range(m):
            for e in range(2):
                c = (1+r)*k[i] + w*e - kp[0]
                J_new[e][i] = sqrt(c) + B*(p[0]*J_old[0][0] + p[1]*J_old[1][0])
                pol[e][:,i] = [k[i], kp[0], c]
            for j in range(m):
                for e in range(2):
                    c = (1+r)*k[i] + w*e - kp[j]
                    if c >= 0: 
                        TMP_J_new = sqrt(c) + B * (p[0]*J_old[0][j] + p[1]*J_old[1][j])
                        if TMP_J_new > J_new[e][i]:
                            J_new[e][i] = TMP_J_new
                            pol[e][:, i] = [k[i], kp[j], c]
        if abs(J_old-J_new).max() < 1e-6:
            break
    return [J_new, pol]

firstcut()
#Gridsize 100: approx 20 sec.