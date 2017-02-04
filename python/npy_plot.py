#! /usr/bin/env python

# plot the colormap and countor for a single data 
import sys
import numpy as np
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt

if __name__ == '__main__':
    data = np.load(sys.argv[1])
    print data.shape
    fig = plt.figure()
    #ax1  = fig.add_subplot(121, aspect='equal')
    ax1  = fig.add_subplot(111, aspect='equal')
    #cax = ax1.imshow(data[1:255,1:255], interpolation='mitchell', origin='lower')
    cax = ax1.imshow(data, interpolation='mitchell', origin='lower')
    ax1.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0, 0.2*data.shape[0], 0.4*data.shape[0], 0.6*data.shape[0], 0.8*data.shape[0], data.shape[0]-1]))
    ax1.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0, 0.2*data.shape[1], 0.4*data.shape[1], 0.6*data.shape[1], 0.8*data.shape[1], data.shape[1]-1]))
    ax1.yaxis.set_ticklabels(['', 0.2, 0.4, 0.6, 0.8, 1])
    ax1.xaxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    for label_i in ax1.get_xticklabels() + ax1.get_yticklabels():
        label_i.set_size(16)
    cb = fig.colorbar(cax)
    for label_i in cb.ax.get_yticklabels():
        label_i.set_size(16)
    #print data[4,:]
    #ax.imshow(data, interpolation='mitchell')
    #ax2  = fig.add_subplot(122, aspect='equal')
    #ax2.contour(data[1:255,1:255])
    #ax2.set_xlim(0, data.shape[1]-1)
    #ax2.set_ylim(0, data.shape[0]-1)

    #ax2 = fig.add_subplot(122)
    #print data[40,:]
    #ax2.plot(data[40,1:255])
    #ax2.set_xlim(0, 255)

    plt.show()
    #fig.set_size_inches(10.0, 10*0.75)
    #plt.savefig("a.pdf")
