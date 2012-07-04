# generates plot of Verification Cases 1 and 2

import numpy as np
import scipy as sp
import matplotlib as ml
import matplotlib.pyplot as plt

# see Table 1
K0    = 0.001       # m s-1;  minimum of hydraulic conductivity
Creep = 1.9014e-4
Cavit = 0.001
A     = 3.1689e-24
g     = 9.81        # m s-2;  acceleration of gravity
n     = 3.0
rhoi  = 910.0       # kg m-3; ice density
rhow  = 1000.0      # kg m-3; freshwater density

c0 = K0 / (rhow * g)

secperday = 3600.0 * 24.0
secpera   = 31556926.0

# see Table 2 for verification case constants

def case1(r):
  alpha  = 3.0e-4
  beta   = 0.0002
  gamma  = 0.15
  H0     = 1000.0
  Wl     = 1.0 + 0.5 * np.cos(alpha * r)
  Wnew   = Wl
  P0     = rhoi * g * H0
  print 'case 1:  constant p_over = %.2f (bar)' % (P0 / 1.0e5)
  pover  = P0 * np.ones(np.size(r))
  P      = P0 * (1.0 - gamma * np.cos(beta * r)**2)
  Y      = 2.0 * beta * r
  sinor  = 2.0 * beta * np.ones(np.size(r))
  sinor[r>0.0] = np.sin(Y[r>0.0]) / r[r>0.0]  # sinor = sin(2 beta r) / r
  divWgradP = Wl * (sinor + 2.0 * beta * np.cos(Y))
  divWgradP -= 0.5 * alpha * np.sin(alpha * r) * np.sin(Y)
  divWgradP *= P0 * gamma * beta
  S = - c0 * rhow * divWgradP
  cvb = Creep * A / Cavit
  vbspeed = cvb * ( pover * gamma * (np.cos(beta*r))**2.0 )**3.0 * Wl
  return Wl, Wnew, pover, P, S, vbspeed

def case2(r):
  alpha0 = 45000.0
  alpha1 = 50000.0
  omega0 = 1.0e-4
  H0 = 2000.0
  n  = 3.0
  Ls = 200000.0
  s  = r / Ls
  ss = s[s<1.0]
  m  = n / (2.0*n + 2.0)
  p  = 1.0 + 1.0/n
  Cs  = H0 / (1.0-1.0/n)**m
  H  = np.zeros(np.size(r))
  H[s<1.0] = Cs * ( p * ss - (1.0/n) + (1.0-ss)**p - ss**p )**m
  Wl   = omega0 * r * np.exp(-r / alpha0)
  Wnew = omega0 * r * np.exp(-r / alpha1)
  Wl[s>=1.0]   = 0.0
  Wnew[s>=1.0] = 0.0
  pover = rhoi * g * H
  P0    = 0.9 * rhoi * g * H0
  print 'case 2:  non-constant p_over with maximum = %.2f (bar)' % (P0 / 1.0e5)
  P       = np.zeros(np.size(r))
  S       = np.zeros(np.size(r))
  vbspeed = np.zeros(np.size(r))
  z = 2.0 * np.pi * r[s<1.0] / (3.0 * Ls)
  y1 = (2.0 - r[s<1.0] / alpha1) * np.sin(z) + z * np.cos(z)
  Cp = 4.0 * np.pi * P0 * omega0 / (9.0 * Ls)
  P[s<1.0] = (P0 / 3.0) * ( 2.0 * np.cos(z) + 1.0 )
  Deltat = 100.0 * secperday
  S[s<1.0] = (rhow / Deltat) * (Wnew[s<1.0] - Wl[s<1.0]) \
             + rhow * c0 * Cp * np.exp(-r[s<1.0] / alpha1) * y1
  y0 = (2.0 - r[s<1.0] / alpha0) * np.sin(z) + z * np.cos(z)
  vbspeed[s<1.0] = - (c0 * Cp / Cavit) * np.exp(-r[s<1.0] / alpha0) * y0 \
                   + (Creep * A / Cavit) * (pover[s<1.0] - P[s<1.0])**n * Wl[s<1.0] \
                   + S[s<1.0] / (rhow * Cavit)
  return Wl, Wnew, pover, P, S, vbspeed

def puttag(ax,tag):
  """
  put a label on a subplot in the lower left corner
  """
  plt.text(0.1, 0.1,tag,weight='bold',fontsize=14.0,
           horizontalalignment='center',verticalalignment='center',
           transform = ax.transAxes)


for case in [1,2]:
  if case==1:
    L = 20000.0
    r = np.arange(0,L * np.sqrt(2.0),50.0)
    Wl, Wnew, pover, P, S, vbspeed = case1(r)
  else:
    L = 250000.0
    r = np.arange(0,L,200.0)
    Wl, Wnew, pover, P, S, vbspeed = case2(r)

  plt.figure(case,figsize=(12.0,7.0)),  plt.clf

  plt.subplot(2,2,1)
  plt.plot(r/1000.0,P/1.0e5,label=r'$P$')
  plt.plot(r/1000.0,pover/1.0e5,'--',label=r'$P_i$')
  plt.grid(True),  plt.legend()
  plt.ylabel('pressure  (bar)')
  puttag(plt.gca(),'(a)')

  plt.subplot(2,2,2)
  if case==1:
    plt.plot(r/1000.0,Wl,label=r'$W_l = W_{l+1}$')
  else:
    plt.plot(r/1000.0,Wnew,label=r'$W_{l+1}$')
    plt.plot(r/1000.0,Wl,'--',label=r'$W_l$')
  plt.grid(True),  plt.legend()
  plt.ylabel('water thickness  (m)')
  puttag(plt.gca(),'(b)')

  plt.subplot(2,2,3)
  plt.plot(r/1000.0,(S / rhow) * secpera,label=r'$S / \rho_w$')
  plt.grid(True),  plt.legend(),  plt.xlabel(r'$r$  (km)')
  plt.ylabel('drainage rate  (m  a-1)')
  puttag(plt.gca(),'(c)')

  plt.subplot(2,2,4)
  plt.plot(r/1000.0,vbspeed * secpera,label=r'$|\mathbf{v}_b|$')
  plt.grid(True),  plt.legend(),  plt.xlabel(r'$r$  (km)')
  plt.ylabel('sliding speed  (m  a-1)')
  puttag(plt.gca(),'(d)')

plt.show()

