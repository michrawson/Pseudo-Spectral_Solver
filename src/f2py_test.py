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

nplots = 24
plotgap = 15

sigma1 = 2.52
sigma2 = 13.4
n = 64+32
h = 15.0/n
dt = 1./2.**4 * 25./32.*13./8.
delta2 = 0.8

x, y, tdata, result = fmain.variable_coeff_wave_eq_pseudo_rk.variable_coeff_wave_eq_pseudo_rk_run(
		nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(h)+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+'-'+format(i, '02')+".png")

sys.exit()
