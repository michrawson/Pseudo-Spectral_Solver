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

if True:
	nplots = 20
	plotgap = 7
	n = 96
	dt = 1./2.**5 * 20./16. * 11./8.
elif False:
	nplots = 20
	plotgap = 8
	n = 64
	dt = 1./2.**5 * 20./16. * 33./32. 
elif False:
	nplots = 20
	plotgap = 4
	n = 160
	dt = 1./2.**5 * 20./16. * 33./32. *.75
elif False:
	nplots = 20
	plotgap = 5
	n = 128
	dt = 1./2.**5 * 20./16. * 33./32. 

sigma1 = 2.*(2.52/2.)**2.
sigma2 = 2.*(13.4/2.)**2.
r = 44.
h = r/n
delta2 = 0.8

x, y, tdata, result = fmain.solve.solve_run(nplots, plotgap, sigma1, sigma2, dt, delta2, h, n)

for i in range(result.shape[0]):
    plt.clf()

    plt.imshow(result[i,:,:], extent=[x[0],x[-1],y[0],y[-1]])

    plt.savefig('wave-'+str(round(r,5))+'-'+str(round(dt,5))+'-'+str(n)+'-'+str(plotgap)+'-'+format(i, '02')+".png")

sys.exit()
