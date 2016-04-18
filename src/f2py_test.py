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

nplots = 16
plotgap = 12

sigma1 = 2.52
sigma2 = 13.4**1.5
n = 64
dt = 1./2.**4
delta2 = 0.8

x, y, tdata, result = fmain.variable_coeff_wave_eq_pseudo_rk.variable_coeff_wave_eq_pseudo_rk_run(
		nplots, plotgap, sigma1, sigma2, dt, delta2, n)

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(round(dt,5))+'-'+str(delta2)+'-'+str(i)+".png")

sys.exit()
