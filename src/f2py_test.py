import numpy as np
import sys
sys.settrace
import fmain
import matplotlib.pyplot as plt

error = np.zeros(9+1,'d')
Nvec = np.asarray((np.ones(9+1,'i')*2) ** range(3,12+1),'i')

# print fmain.__doc__
# print
# print fmain.finite_differences.__doc__
# print
# print fmain.variable_coeff_wave_eq.__doc__

# fmain.finite_differences.finite_differences_run(error)

# print 'error',error


# x, tdata, result = fmain.variable_coeff_wave_eq.variable_coeff_wave_eq_run()
x, tdata, result = fmain.variable_coeff_wave_eq_pseudo.variable_coeff_wave_eq_pseudo_run()
# print 'x',x.shape
# print 'x',x
# print 'tdata',tdata.shape
# print 'tdata',tdata
# print 'result',result.shape
# print 'result',result[0,:]
# print 'result',result[1,:]
# print 'result',result[2,:]

for i in range(result.shape[0]):
    plt.clf() 
    ax = plt.figure().add_subplot(1,1,1)
    ax.plot(x, result[i,:])
    plt.ylim(-.1,1)
    plt.savefig('wave'+str(i)+".png")

sys.exit()
