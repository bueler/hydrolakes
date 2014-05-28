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

def criticalW(Po,vb):
  sbcube = c1 * vb / (c2 * A)
  Wc = Wr * (sbcube / (sbcube + Po**3.0))
  return Wc

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

fig = plt.figure(1,figsize=(6.0,4.0))
H = 1000.0
Po = rhoi * g * H
vb = [0.0, 10.0, 100.0, 1000.0]
for j in range(4):
  P = psteady(W,Po,vb[j]/spera)
  Wc = criticalW(Po,vb[j]/spera)
  plt.plot(Wc,0.0,'ko',markersize=8.0,markeredgecolor='k',markerfacecolor='k')
  plt.hold(True)
  plt.plot(W[W>=0.95*Wc],P[W>=0.95*Wc]/1.0e5,'k',lw=2.0)
  plt.text(W[40 + 40*j],(P[40 + 40*j]/1.0e5) - 7.0 + 5.0 * j,r'$|\mathbf{v}_b| =$ %d m/a' % vb[j],rotation=8.0*j)
  if j == 1:
    plt.plot([0.0, W[1]],[0.0, P[1]/1.0e5],'k',lw=2.0)
plt.plot(W,Po * (W/Wr)**3.5 / 1.0e5,'k--',lw=2.5)
plt.hold(False)

#plt.gca().set_aspect('equal')
#plt.gca().autoscale(tight=True)
plt.gca().set_clip_on(False)
plt.xlabel('W  (m)')
plt.ylabel('P  (bar)')
plt.axis([-0.02, 1.19, -3.0, 110.0])
plt.tight_layout()

#plt.show()
fig.savefig('psteady-vb.pdf', bbox_inches='tight')
