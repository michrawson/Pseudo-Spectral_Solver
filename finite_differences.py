import numpy as np
import finite_differences
import matplotlib.pyplot as plt

error = np.zeros((2**(9+2)),'d')

print finite_differences.__doc__
print finite_differences.run.__doc__

finite_differences.run(error)

print error

plt.plot(error[0:20])
plt.show()
