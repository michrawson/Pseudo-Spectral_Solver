import numpy as np
import sys
import finite_differences
import matplotlib.pyplot as plt

error = np.zeros(9,'d')
Nvec = (np.ones(9+1,'d')*2) ** range(3,12+1)

print finite_differences.__doc__
print finite_differences.run.__doc__

finite_differences.run(error)

print error
print 
print Nvec

ax = plt.figure().add_subplot(1,1,1)
ax.plot(error[0:20])
ax.plot(Nvec)
ax.set_yscale('log')
plt.show()
