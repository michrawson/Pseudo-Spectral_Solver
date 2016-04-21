import numpy as np
import sys
sys.settrace
import fmain
import matplotlib.pyplot as plt

# print fmain.__doc__
# print
# print fmain.finite_differences.__doc__
# print
# print fmain.variable_coeff_wave_eq.__doc__

nplots = 20
plotgap = 16

sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
n = 128
r = 50.
h = r/n
dt = 1./2.**5 * 20./16.
delta2 = 0.8

x, y, tdata, result = fmain.variable_coeff_wave_eq_pseudo_rk.variable_coeff_wave_eq_pseudo_rk_run(
		nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+'-'+format(i, '02')+".png")

sys.exit()
