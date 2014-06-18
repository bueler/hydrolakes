#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

spera = 31556926.0

e0    = 0.69          # reference void ratio at N0; pure
N0    = 1.0e3         # reference effective pressure; Pa
Cc    = 0.12          # compressibility; pure

g     = 9.81          # m s-2
rhoi  = 910.0         # kg m-3

H     = 1000.0        # m; typical thickness

Po = rhoi * g * H

N = np.linspace(0.1*N0,1.5*Po,500)
e = e0 - Cc * np.log10(N / N0)   # eqn (3b) of Tylaczyk et al (2000b)

# FIGURE (a): Original Tulaczyk
fig = plt.figure(1,figsize=(5.0,7.0))
plt.subplot(2,1,1)
plt.plot(e,N / 1.0e5,'k',lw=2.5)
plt.plot(e0,N0 / 1.0e5,'ko',markersize=8.0,markeredgecolor='k',markerfacecolor='k')
plt.text(0.9*e0,600.0*N0 / 1.0e5,r'$(e_0,N_0)$',fontsize=14.0)
plt.xlabel(r'$e$',fontsize=16.0)
plt.ylabel(r'$N_{til}$  (bar)',fontsize=14.0)
plt.axis([0.0, 0.8, -2.5, 130.0])
plt.text(-0.13,1.35*Po / 1.0e5,'(a)',fontsize=16.0)
plt.tight_layout()


# FIGURE (b): Redo with:
#   (1) e replaced by Wtil:
#           e = V_w / V_s = (W_til dx dy) / (eta dx dy) = W_til / eta
#   where eta is fully-compressed thickness of solids
#   (2) but  W_til^max = e^max * eta  so eta can be determined from W_til^max
#   and e^max
#   (3) and  N_min = delta P_o  and the relation
#       e^max = e0 - Cc * log10(N_min / N0)
#   gives
#                        W_til^max
#       eta = ------------------------------------
#              e0 - Cc * log10(delta P_o / N0)
#   (4) finally N_til <= P_o
Wtilmax = 2.0   # m
delta   = 0.02
eta     = Wtilmax / (e0 - Cc * np.log10(delta * Po / N0))

plt.subplot(2,1,2)
# background dotted (and smaller)
plt.plot(e * eta,N / 1.0e5,'k:',lw=2.0)
# new solid part
N = np.linspace(delta*Po,Po,500)
e = e0 - Cc * np.log10(N / N0)
Wtil = e * eta
plt.plot(Wtil,N / 1.0e5,'k',lw=2.5)
Wtilflat = np.array([0.0,min(Wtil)])
Poflat = Po * np.ones(np.shape(Wtilflat))
plt.plot(Wtilflat,Poflat / 1.0e5,'k',lw=2.5)
# annotations
plt.text(0.2 * np.min(Wtil),1.05*Po / 1.0e5,r'$N_{til} = P_o$',fontsize=14.0)
plt.plot(Wtilmax,delta*Po / 1.0e5,'ko',markersize=8.0,markeredgecolor='k',markerfacecolor='k')
plt.text(0.9*Wtilmax,5.0*delta*Po / 1.0e5,r'$(W_{til}^{max},\delta P_o)$',fontsize=14.0)
plt.xlabel(r'$W_{til}$  (m)',fontsize=14.0)
plt.ylabel(r'$N_{til}$  (bar)',fontsize=14.0)
plt.axis([0.0*eta, 0.8*eta, -2.5, 130.0])
plt.text(-0.13*eta,1.35*Po / 1.0e5,'(b)',fontsize=16.0)
plt.tight_layout()
#plt.show()
fig.savefig('Ntilfunctions.pdf', bbox_inches='tight')


