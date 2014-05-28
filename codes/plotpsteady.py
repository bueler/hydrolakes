#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

spera = 31556926.0

A     = 3.1689e-24    # ice softness (Pa-3 s-1)
c1    = 0.500         # m-1
c2    = 0.040         # [pure]
g     = 9.81          # m s-2
rhoi  = 910.0         # kg m-3
Wr    = 1.0           # m

Y0    = 0.001         # regularization; m

def psteady(W,Po,vb):
  '''Computes P(W) in steady state.'''
  if np.any(Po < 0.0):
    print 'psteady() requires nonnegative overburden pressure Po'
    sys.exit(1)
  if np.any(vb < 0.0):
    print 'psteady() requires nonnegative sliding speed vb'
    sys.exit(2)
  sbcube = c1 * vb / (c2 * A)
  frac = np.maximum(0.0, Wr - W) / (W + Y0)
  P = Po - (sbcube * frac)**(1.0/3.0)
  P[P < 0.0] = 0.0
  return P

W = np.linspace(0.0,1.2*Wr,501)

plt.figure(1)
H = 1000.0
Po = rhoi * g * H
for vb in np.array([0.0, 10.0, 100.0, 1000.0]) / spera:
  P = psteady(W,Po,vb)
  plt.plot(W,P/1.0e5,'k',lw=2.0)
  plt.hold(True)
  plt.plot(min(W[P>0.0]),0.0,'ko',markersize=12.0,markeredgecolor='k',markerfacecolor='k')
plt.hold(False)

#gca().set_aspect('equal')
plt.gca().autoscale(tight=True)
plt.xlabel('W  (m)')
plt.ylabel('P  (bar)')
plt.axis([0.0, 1.19, 0.0, 110.0])

plt.show()

