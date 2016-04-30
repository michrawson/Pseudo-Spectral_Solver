import numpy as np
import sys
sys.settrace
import fmain
import matplotlib.pyplot as plt

#print fmain.__doc__
#print
#print fmain.solve.__doc__
#print
#print fmain.solve.solve_run.__doc__

# nplots = 20
# plotgap = 36
# n = 64
# dt = 1./2.**4 * 23./16.

nplots = 20
plotgap = 21
n = 96
dt = 1./2.**5 * 20./16. * 11./8.

# nplots = 20
# plotgap = 14
# n = 128
# dt = 1./2.**5 * 20./16.

# nplots = 20
# plotgap = 10
# n = 160
# dt = 1./2.**5 * 20./16. * 33./32. *.75


sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
r = 44.
h = r/n
delta2 = 0.8

x, y, tdata, result = fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+'-'+format(i, '02')+".png")

sys.exit()
