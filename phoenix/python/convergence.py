# generates convergence plot for wpi.c on Verification Case 1

"""
bueler@bueler-lemur ~/icerepo/UAF-misc/betterhydro/petsc$ svnversion 
1706M
bueler@bueler-lemur ~/icerepo/UAF-misc/betterhydro/petsc$ for M in 21 41 81 161 321 641; do     echo "case M=$M:";     mpiexec -n 4 ./wpi -fd -wpi_case 1 -da_grid_x $M -da_grid_y $M |grep errors;   done
case M=21:
6:errors:    av |psi-exact| = 1.411e+04 (Pa),   max |psi-exact| / max |exact| = 7.789e-03
9:errors:    av |Wnew-exact| = 4.032e-03 (m),   max |Wnew-exact| / max |exact| = 1.863e-02
case M=41:
6:errors:    av |psi-exact| = 3.926e+03 (Pa),   max |psi-exact| / max |exact| = 2.153e-03
9:errors:    av |Wnew-exact| = 2.667e-03 (m),   max |Wnew-exact| / max |exact| = 1.006e-02
case M=81:
6:errors:    av |psi-exact| = 9.914e+02 (Pa),   max |psi-exact| / max |exact| = 5.304e-04
9:errors:    av |Wnew-exact| = 1.504e-03 (m),   max |Wnew-exact| / max |exact| = 5.536e-03
case M=161:
6:errors:    av |psi-exact| = 2.416e+02 (Pa),   max |psi-exact| / max |exact| = 1.265e-04
9:errors:    av |Wnew-exact| = 7.990e-04 (m),   max |Wnew-exact| / max |exact| = 2.872e-03
case M=321:
6:errors:    av |psi-exact| = 6.368e+01 (Pa),   max |psi-exact| / max |exact| = 3.316e-05
9:errors:    av |Wnew-exact| = 4.092e-04 (m),   max |Wnew-exact| / max |exact| = 1.403e-03
case M=641:
6:errors:    av |psi-exact| = 1.535e+01 (Pa),   max |psi-exact| / max |exact| = 7.957e-06
9:errors:    av |Wnew-exact| = 2.079e-04 (m),   max |Wnew-exact| / max |exact| = 2.555e-03
"""

import numpy as np
import scipy as sp
import matplotlib as ml
import matplotlib.pyplot as plt

err = np.array([[ 21.0, 1.411e+04, 7.789e-03, 4.032e-03, 1.863e-02],
                [ 41.0, 3.926e+03, 2.153e-03, 2.667e-03, 1.006e-02],
                [ 81.0, 9.914e+02, 5.304e-04, 1.504e-03, 5.536e-03],
                [161.0, 2.416e+02, 1.265e-04, 7.990e-04, 2.872e-03],
                [321.0, 6.368e+01, 3.316e-05, 4.092e-04, 1.403e-03],
                [641.0, 1.535e+01, 7.957e-06, 2.079e-04, 2.555e-03]])

dx = 40000.0 / (err[:,0]-1)
print dx

plt.figure(1,figsize=(10.0,6.0)),  plt.clf

plt.subplot(2,1,1)
p = np.polyfit(np.log(dx),np.log(err[:,1]),1)
print 'rate for av  psi  = dx^%.2f' % p[0]
plt.loglog(dx, err[:,1],'bo',markersize=10,
           label=r'av $\psi$ error: $O(\Delta x^{%.2f})$' % p[0])
plt.loglog(dx, np.exp(np.polyval(p,np.log(dx))), 'r--')
maxerr = err[:,2] * 89.271e5  # put back in Pa
p = np.polyfit(np.log(dx),np.log(maxerr),1)
print 'rate for max psi  = dx^%.2f' % p[0]
plt.loglog(dx, maxerr,'gs',markersize=10,
           label=r'max $\psi$ error: $O(\Delta x^{%.2f})$' % p[0])
plt.loglog(dx, np.exp(np.polyval(p,np.log(dx))), 'r:')
plt.grid(True)
plt.ylabel('error  (Pa)')
plt.legend(loc='upper left')
ax = plt.gca()
ax.set_xticks(dx)
#ax.set_xticklabels(('2 km', '1 km', '500 m', '250 m', '125 m', '62.5 m'))
ax.set_xticklabels((' ',' ',' ',' ',' ',' '))
plt.axis([40.0, 3000.0, 1.0e1, 2.0e5])

plt.subplot(2,1,2)
p = np.polyfit(np.log(dx),np.log(err[:,3]),1)
print 'rate for av  Wnew = dx^%.2f' % p[0]
plt.loglog(dx, err[:,3],'bo',markersize=10,
           label=r'av $W^{l+1}$ error: $O(\Delta x^{%.2f})$' % p[0])
plt.loglog(dx, np.exp(np.polyval(p,np.log(dx))), 'r--')
p = np.polyfit(np.log(dx[0:5]),np.log(err[0:5,4]),1)
print 'rate for max Wnew = dx^%.2f' % p[0]
plt.loglog(dx, err[:,4],'gs',markersize=10,
           label=r'max $W^{l+1}$ error: $O(\Delta x^{%.2f})$' % p[0])
plt.loglog(dx, np.exp(np.polyval(p,np.log(dx))), 'r:')
plt.grid(True)
plt.xlabel(r'$\Delta x = \Delta y$  (m)')
ax = plt.gca()
ax.set_xticks(dx)
ax.set_xticklabels(('2 km', '1 km', '500 m', '250 m', '125 m', '62.5 m'))
plt.axis([40.0, 3000.0, 1.0e-4, 1.0e-1])
plt.ylabel('error  (m)')
plt.legend(loc='upper left')

plt.show()

