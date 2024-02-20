# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as si

def duffing (y,t, g, a, b, A, omega):
 return np.array([y[1], -2*g*y[1]-(a+b*y[0]*y[0])*y[0]\
   +A*np.cos(omega*t)])

def solve (A, ttrans, tfinal):
  y0 = np.zeros(2)
  h = 1e-3
  t = np.arange (0, tfinal+h, h)
  sol = si.odeint (duffing, y0, t, args=(5e-2, -1, 1, A, 1.4))
  i = t>=ttrans
  return (t[i], sol[i])

def section (A, t0, N):
  y0 = np.zeros(2)
  t = np.array([0., t0])
  sol0 = si.odeint (duffing, y0, t, args=(5e-2, -1, 1, A, 1.4), h0=0.001,\
    mxstep=int(2000*t0))
  T = 2*np.pi/1.4
  t = t0+np.arange(N)*T
  sol = si.odeint (duffing, sol0[-1], t, args=(5e-2, -1, 1, A, 1.4),\
    h0=0.001, mxstep=int(2000*(N-1)*T))
  return (sol)

# Q5
# Cas 1: solution périodique avec période T
# Cas 2: solution périodique avec période 2T
# Cas 3: solution périodique avec période 4T
# Cas 4: solution chaotique

# Q8
# Cas 1: Un seul point (aux erreurs numériques près) : la
#        solution revient au même point après une période
# Cas 2: Deux points (aux erreurs numériques près) : la
#        solution oscille entre deux valeurs différentes à chaque
#        période. La période de l'oscillateur est doublée.
# Cas 3: Quatre points (aux erreurs numériques près) : la
#        solution oscille entre quatre valeurs différentes à chaque
#        période. La période de l'oscillateur est quadruplée.
# Cas 4: Chaque point est différent, la solution est chaotique.

A=np.array ([0.1,0.321,0.338,0.35])
ttrans=np.array([150.,750.,1800.,2800])
tfinal=np.array([200.,800.,2000.,5000.])
for i in range(4):
  sol = solve (A[i], ttrans[i], tfinal[i])
  plt.plot (sol[0], sol[1][:,0], label="x(t)")
  plt.plot (sol[0], sol[1][:,1], label="v(t)")
  plt.title("A = " + str(A[i]))
  plt.xlabel("t")
  plt.legend ()
  plt.savefig("xvt-" + str(i+1) + ".png")
  plt.show ()
  plt.plot (sol[1][:,0], sol[1][:,1])
  plt.title("A = " + str(A[i]))
  plt.xlabel("x")
  plt.ylabel("v")
  plt.savefig("phase-" + str(i+1) + ".png")
  plt.show()
  sec = section (A[i], tfinal[i], 20)
  plt.plot (sec[:,0], sec[:,1], "o")
  plt.title("Section de Poincaré pour A = " + str(A[i]))
  plt.xlabel("x")
  plt.ylabel("v")
  plt.xlim (-2,2)
  plt.ylim (-1.5,1.5)
  plt.savefig("section-" + str(i+1) + ".png")
  plt.show()
