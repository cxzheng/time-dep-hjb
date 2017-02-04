#! /usr/bin/env python

# =============================================================================
# Given the ground truth solution (high res.) and the low resolution solution,
# it reports the L_1 and L_{\infty}
# NOTICE: the high resolution should be multiple of the low resolution
# =============================================================================

import sys
import numpy as np
from termcolor import colored
#import matplotlib.pyplot as plt

def die(msg):
    sys.stderr.write(colored('ERROR: ', 'red', attrs=['bold']))
    print >> sys.stderr, colored(msg, 'yellow')
    sys.exit(1)

# plot the data
def plot(data):
    fig = plt.figure()
    ax1 = fig.add_subplot(121, aspect='equal')
    cax = ax1.imshow(data, interpolation='bilinear', origin='lower') #, vmax=0.0094)
    fig.colorbar(cax)
    ax2  = fig.add_subplot(122, aspect='equal')
    ax2.contour(data)
    ax2.set_xlim(0, data.shape[1])
    ax2.set_ylim(0, data.shape[0])
    plt.show()

def find_max(data):
    maxv = 0
    for iy in xrange(data.shape[0]):
        for ix in xrange(data.shape[1]):
            if data[iy,ix] > maxv:
                maxv = data[iy,ix]
                coord = (ix, iy)
    print "Val = %f" % maxv, "  Max coord: %d, %d" % coord

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print >> sys.stderr, "Usage: %s <ground truth data> <approx data>" % sys.argv[0]
        die("Invalid arguments")

    # load data and check their dimensions
    gData = np.load(sys.argv[1])
    if len(gData.shape) != 2:
        die("%s: 2D data is expected" % sys.argv[1])
    cData = np.load(sys.argv[2])
    if len(cData.shape) != 2:
        die("%s: 2D data is expected" % sys.argv[2])
    if (gData.shape[0]-1) % (cData.shape[0]-1) or (gData.shape[1]-1) % (cData.shape[1]-1):
        die("high res. is not a multiple of low res")

    # create low res very of the array   
    dy = (gData.shape[0]-1)/(cData.shape[0]-1)
    dx = (gData.shape[1]-1)/(cData.shape[1]-1)
    #print dx, dy
    iData = np.empty(cData.shape)
    for iy in xrange(cData.shape[0]):
        for ix in xrange(cData.shape[1]):
            iData[iy,ix] = gData[iy*dy, ix*dx]
    #plot(iData)

    # compute L_0 and L_{\infty} error
    err = np.abs(iData - cData)
    #find_max(err)
    #print "error:", np.mean(err), np.max(err)
    print np.mean(err), np.max(err)
    #plot(err)
