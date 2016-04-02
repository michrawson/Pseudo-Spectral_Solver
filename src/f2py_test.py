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


p_xi = fmain.zero_finder.triginterp_caller(np.pi)
print p_xi

sys.exit()

x, y, tdata, result = fmain.variable_coeff_wave_eq_pseudo_rk.variable_coeff_wave_eq_pseudo_rk_run()

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave'+str(i)+".png")

sys.exit()
