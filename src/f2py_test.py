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

fmain.finite_differences.finite_differences_run(error)

print 'error',error


x, tdata, result = fmain.variable_coeff_wave_eq.variable_coeff_wave_eq_run()
print 'x',x.shape
print 'x',x
print 'tdata',tdata.shape
print 'tdata',tdata
print 'result',result.shape
# print 'result',result

for i in range(result.shape[0]):
    print result[i,:]

# print result[0,:].T - result[20,:].T

ax = plt.figure().add_subplot(1,1,1)
# ax.plot(x,result[0,:].T, 'r--',
#         x, result[1,:].T, 'bs',
#         x, result[2,:].T, 'g^')
        # x, result[3,:].T, 'y--',
        # x, result[4,:].T, 'os',
        # x, result[5,:].T, 'b^',
        # x, result[6,:].T, 'p--')

ax.plot(x,result[0,:])
# ax.plot(x,result[20,:].T)
# ax.plot(x,result[30,:].T)
# ax.plot(x,result[40,:].T)
# ax.plot(x,result[50,:].T)
# ax.loglog(Nvec, error)
# ax.loglog(Nvec, Nvec**(-4.),'--')
plt.show()
