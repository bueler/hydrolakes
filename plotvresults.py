#!/usr/bin/env python

from pylab import *
from sys import exit
import numpy as np

# PLOTVRESULTS  semilog plot to show how v depends on grid spacing

# data from 'results.txt'
dxkm = np.array([100., 50., 25., 15., 10., 5.])
vmaxpera = np.array([694.514, 2240.086, 4378.39, 6887.370, 9992.122, 12847.341])

semilogy(dxkm, vmaxpera, 'o', markersize=10)

axis([0, 105, 500, 20000])
ax = subplot(111)

xlabel(r'$\Delta x$  (km)', fontsize=16)
ax.set_xticks(dxkm.tolist() + [0.])
ax.set_xticklabels(('100','50','25','15','10','5','0'))

ylabel(r'$\max |\mathbf{v}|$  (m a-1)', fontsize=16)
ax.set_yticks([1000, 2000, 5000, 10000])
ax.set_yticklabels(('1000','2000','5000','10000'))

grid('on', linewidth=2.0)

show()
