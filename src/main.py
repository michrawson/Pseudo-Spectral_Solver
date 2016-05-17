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

# MASTER
nplots = 25
plotgap = 5
n = 96
dt = 1./2.**5 * 20./16. * 11./8.
tmax = nplots*plotgap*dt

# RK5
# nplots = 20
# plotgap = 14
# n = 96
# dt = 1./2.**5 * 20./16. * 11./8. * 5./4.


# RK4

# nplots = 20
# plotgap = 9
n = 64
dt = 1./2.**4 * 23./16.
plotgap = tmax/dt/nplots

sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
r = 50.
h = r/n
delta2 = 0.8

x, y, tdata, result, int_res, int2_res, max_res, min_res = \
    fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

# print 'max_res',max_res

plt.clf()
plt.plot(int_res)
plt.title('integral u(n) '+str(n)+'x'+str(n))
plt.savefig('int_res-'+str(n)+".png")

plt.clf()
plt.plot(int2_res)
plt.title('integral u(n)^2 ' +str(n)+'x'+str(n))
plt.savefig('int2_res-'+str(n)+".png")

plt.clf()
plt.plot(max_res)
plt.title('max u(n) '+str(n)+'x'+str(n))
plt.savefig('max_res-'+str(n)+".png")

plt.clf()
plt.plot(min_res)
plt.title('min u(n) '+str(n)+'x'+str(n))
plt.savefig('min_res-'+str(n)+".png")

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], vmin=-.1, vmax=1.1, extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(n)+'-'+format(i, '02')+".png")


sys.exit()
###################################################################


# nplots = 20
plotgap = 5
n = 96
dt = 1./2.**5 * 20./16. * 11./8.

sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
r = 50.
h = r/n
delta2 = 0.8

x, y, tdata, result, int_res, int2_res, max_res, min_res = \
    fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

plt.clf()
plt.plot(int_res)
plt.title('integral u(n) '+str(n)+'x'+str(n))
plt.savefig('int_res-'+str(n)+".png")

plt.clf()
plt.plot(int2_res)
plt.title('integral u(n)^2 ' +str(n)+'x'+str(n))
plt.savefig('int2_res-'+str(n)+".png")

plt.clf()
plt.plot(max_res)
plt.title('max u(n) '+str(n)+'x'+str(n))
plt.savefig('max_res-'+str(n)+".png")

plt.clf()
plt.plot(min_res)
plt.title('min u(n) '+str(n)+'x'+str(n))
plt.savefig('min_res-'+str(n)+".png")

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], vmin=-.1, vmax=1.1, extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(n)+'-'+format(i, '02')+".png")


###################################################################


nplots = 10
# plotgap = 3
n = 128
dt = 1./2.**5 * 20./16.
plotgap = tmax/dt/nplots

sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
r = 50.
h = r/n
delta2 = 0.8

x, y, tdata, result, int_res, int2_res, max_res, min_res = \
    fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

plt.clf()
plt.plot(int_res)
plt.title('integral u(n) '+str(n)+'x'+str(n))
plt.savefig('int_res-'+str(n)+".png")

plt.clf()
plt.plot(int2_res)
plt.title('integral u(n)^2 ' +str(n)+'x'+str(n))
plt.savefig('int2_res-'+str(n)+".png")

plt.clf()
plt.plot(max_res)
plt.title('max u(n) '+str(n)+'x'+str(n))
plt.savefig('max_res-'+str(n)+".png")

plt.clf()
plt.plot(min_res)
plt.title('min u(n) '+str(n)+'x'+str(n))
plt.savefig('min_res-'+str(n)+".png")

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], vmin=-.1, vmax=1.1, extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(n)+'-'+format(i, '02')+".png")

sys.exit()
###################################################################

# nplots = 80
# plotgap = 3
n = 160
dt = 1./2.**5 * 20./16. * 33./32. *.75
plotgap = tmax/dt/nplots

sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
r = 50.
h = r/n
delta2 = 0.8

x, y, tdata, result, int_res, int2_res, max_res, min_res = \
    fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

plt.clf()
plt.plot(int_res)
plt.title('integral u(n) '+str(n)+'x'+str(n))
plt.savefig('int_res-'+str(n)+".png")

plt.clf()
plt.plot(int2_res)
plt.title('integral u(n)^2 ' +str(n)+'x'+str(n))
plt.savefig('int2_res-'+str(n)+".png")

plt.clf()
plt.plot(max_res)
plt.title('max u(n) '+str(n)+'x'+str(n))
plt.savefig('max_res-'+str(n)+".png")

plt.clf()
plt.plot(min_res)
plt.title('min u(n) '+str(n)+'x'+str(n))
plt.savefig('min_res-'+str(n)+".png")


# for i in range(result.shape[0]):
#     plt.clf()
#
#     plt.imshow(result[i,:,:], vmin=-.1, vmax=1.1, extent=[x[0],x[-1],y[0],y[-1]])
#
#     plt.savefig('wave-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+'-'+format(i, '02')+".png")


sys.exit()
