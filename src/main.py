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

# RK5
# nplots = 20
# plotgap = 14
# n = 96
# dt = 1./2.**5 * 20./16. * 11./8. * 5./4.


# RK4

# nplots = 20
# plotgap = 36
# n = 64
# dt = 1./2.**4 * 23./16.

nplots = 40
plotgap = 10
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
r = 50.
h = r/n
delta2 = 0.8

x, y, tdata, result, int_res, int2_res, max_res, min_res = \
    fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], vmin=-.1, vmax=1.1, extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+'-'+format(i, '02')+".png")

plt.clf()
plt.plot(int_res)
plt.title('integral u(n)')
plt.savefig('int_res-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+".png")

plt.clf()
plt.plot(int2_res)
plt.title('integral u(n)^2')
plt.savefig('int2_res-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+".png")

plt.clf()
plt.plot(max_res)
plt.title('max u(n)')
plt.savefig('max_res-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+".png")

plt.clf()
plt.plot(min_res)
plt.title('min u(n)')
plt.savefig('min_res-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+".png")

sys.exit()
