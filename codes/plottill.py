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

N = np.linspace(0.1*N0,1.2*Po,500)
Ntall = np.linspace(1.2*Po,1.5*Po,100) 

e = e0 - Cc * np.log10(N / N0)   # eqn (3b) of Tylaczyk et al (2000b)

etall = e0 - Cc * np.log10(Ntall / N0)

eflat = np.array([0.01, 0.79])
Poflat = Po * np.ones(np.shape(eflat))

fig = plt.figure(1,figsize=(5.0,4.0))

plt.plot(e,N / 1.0e5,'k',lw=2.0)
plt.plot(etall,Ntall / 1.0e5,'k:',lw=2.0)

plt.plot(eflat,Poflat / 1.0e5,'k--',lw=2.0)
plt.text(e0,1.05*Po / 1.0e5,r'$P_o$',fontsize=14.0)

plt.plot(e0,N0 / 1.0e5,'ko',markersize=6.0,markeredgecolor='k',markerfacecolor='k')
plt.text(0.9*e0,600.0*N0 / 1.0e5,r'$(e_0,N_0)$',fontsize=14.0)

plt.hold(False)
plt.xlabel(r'$e$')
plt.ylabel(r'$N_{til}$  (bar)')
plt.axis([0.0, 0.8, -2.5, 130.0])
plt.tight_layout()
plt.show()
#fig.savefig('foo.pdf', bbox_inches='tight')


