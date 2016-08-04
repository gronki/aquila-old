#!/bin/python

import numpy as np
import matplotlib
matplotlib.use('GTK3Cairo')
import matplotlib.pyplot as plt

d = np.loadtxt('temp/aa',skiprows=2)
plt.plot(d[:,1],d[:,2],'k+')
plt.show()
