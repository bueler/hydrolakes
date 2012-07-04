import numpy as np
import scipy as sp
import matplotlib as ml
import matplotlib.pyplot as plt
import scipy.integrate as spi

rhoi  = 910.0       # kg m-3; ice density
rhow  = 1000.0      # kg m-3; freshwater density
g     = 9.81        # m s-2;  acceleration of gravity
sigma = 7.0/2.0     # pure; from Flowers & Clarke (2002)

Kmax  = 0.01        # m s-1;  maximum of hydraulic conductivity
Kmin  = 0.0001      # m s-1;  minimum of hydraulic conductivity
K0    = sp.sqrt(Kmax * Kmin) # geometric average

Y0    = 1.0   # m
S0    = rhow * (1.0/31556926.0) # = 1 m a-1 added to water layer, in kg m-2 s-1

H0    = 1000.0
L     = 20000.0

profileC = H0 / (2.0**(3.0/8.0))

ell = L - 1000.0


def phi(s):
  """inside profile function"""
  return 4.0 * s - 3.0 * s**(4.0/3.0) + 3.0 * (1.0 - s)**(4.0/3.0) - 1.0

def HH(x):
  """compute the 'Bueler profile', Greve & Blatter (2009) subsection 5.6.3"""
  s = abs(x / L)  # scaled coordinate
  if s <= 1.0:
    return profileC * phi(s)**(3.0/8.0)
  else:
    return 0.0

def dHHdx(x):
  """gradient of the profile"""
  s = abs(x / L)  # scaled coordinate
  if s < 1.0:
    CC = (3.0 * profileC) / (2.0 * L) 
    return CC * phi(s)**(-5.0/8.0) * (1.0 - s**(1.0/3.0) + (1.0 - s)**(1.0/3.0))
  else:
    return 0.0

def pi(x):
  """overburden pressure"""
  return rhoi * g * HH(x)

def dpidx(x):
  """gradient of overburden pressure"""
  return rhoi * g * dHHdx(x)

# constant in ODE
Omega = S0 * g * (Y0**sigma) / (3.0 * ell**2.0 * K0)

def func(W, x):
  """right-hand-side of ODE"""
  Wsig = W**sigma
  return (- Omega * x**3.0 - dpidx(x) * Wsig * W) / (sigma * pi(x) * Wsig)

# setup grid
N = 190
dx = ell / N
x = np.arange(0.0,ell+dx,dx)  # 0 <= x <= ell
xstag = x[:-2] + (dx/2.0)  # needed to show Q

# solve ODE for W(x) on grid, with decreasing x because we know W(x=ell)
W_at_ell = Y0
#W = spi.odeint(func, W_at_ell, x[::-1])
#W = W.flatten()
#W = W[::-1]
#print '(x,W) near x=l:\n', x[-3:], '\n', W[-3:]
#print 'max W = %f (m)' % max(W)

# we believe: P = rhow * g * H * W**sigma / Y0**sigma
# solve ODE for P(x), with decreasing x because we know P(x=ell)
POmega = S0 * g / (3.0 * ell**2.0 * K0 * Y0)
def Pfunc(P, x):
  """right-hand-side of ODE for P(x)"""
  frac = (pi(x) / P)**(1.0/sigma)
  return - POmega * x**3.0 * frac
P_at_ell = pi(ell) * W_at_ell**sigma / (Y0**sigma)
P = spi.odeint(Pfunc, P_at_ell, x[::-1])
P = P.flatten()
P = P[::-1]
print 'max P = %f (bar)' % max(P/1.0e5)

# compute H for plot (which goes past ell)
H = np.zeros(len(x))
for j in range(len(x)):
  H[j] = HH(x[j])
xextend = np.arange(ell,L+10.0,10.0)  # ell <= x <= L
Hextend = np.zeros(len(xextend))
for j in range(len(xextend)):
  Hextend[j] = HH(xextend[j])

# compute W
W = np.zeros(len(x))
for j in range(len(x)):
  pij = pi(x[j])
  W[j] = Y0 * (pij / P[j])**(1.0/sigma) 

# compute S for plot
S = S0 * (x / ell)**2.0

plt.figure(1)

plt.subplot(311)
plt.plot(x/1000.0,H,xextend/1000.0,Hextend,'--')
plt.axis( [0.0, L/1000.0, 0.0, 1.1 * H[0]] ) # [xmin, xmax, ymin, ymax]
plt.ylabel('H  (m)')
plt.title('inputs')

plt.subplot(312)
plt.plot(x/1000.0,rhoi * g * H / 1.0e5)
plt.axis( [0.0, L/1000.0, 0.0, 1.1 * rhoi * g * H[0] / 1.0e5] )
plt.ylabel('overburden  (bar)')

plt.subplot(313)
plt.plot(x/1000.0,S / 31556926.0)
plt.ylabel('S  (kg m-2 a-1)')
plt.xlabel('x  (km)')

plt.figure(2)

plt.subplot(311)
plt.plot(x/1000.0,W)
plt.axis( [0.0, L/1000.0, 0.0, 1.1 * W[0]] )
plt.ylabel('W  (m)')
plt.title('outputs')

Wstag = (W[:-2] + W[1:-1]) / 2.0
Q = - (K0 / (rhow * g)) * Wstag * (P[1:-1] - P[:-2]) / dx  # staggered vals
Qexact = (S0 / (3.0 * rhow * ell**2.0)) * xstag**3
plt.subplot(312)
plt.plot(xstag/1000.0,Q,xstag/1000.0,Qexact,'--')
plt.axis( [0.0, L/1000.0, 0.0, 1.1 * max(Q)] )
plt.ylabel('Q  (m2 s-1)')

plt.subplot(313)
plt.plot(x/1000.0,P/1.0e5)
plt.axis( [0.0, L/1000.0, 0.0, 1.1 * max(P)/1.0e5] )
plt.ylabel('P  (bar)')
plt.xlabel('x  (km)')

plt.show()

