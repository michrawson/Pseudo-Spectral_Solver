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

x, tdata, result = fmain.variable_coeff_wave_eq_pseudo_rk.variable_coeff_wave_eq_pseudo_rk_run()

for i in range(result.shape[0]):
    plt.clf() 
    ax = plt.figure().add_subplot(1,1,1)
    ax.plot(x, result[i,:])
    plt.ylim(-.1,1)
    plt.savefig('wave'+str(i)+".png")

sys.exit()
