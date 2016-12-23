#!/usr/bin/python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#           Name: ksTest.py
#         Author: Ra Inta, 20140111
#  Last Modified: 20160519
#    Description: This takes a search_histogram.txt.$jobNumber
# and performs a Kolmogorov-Smirnoff test by comparing the cumulative distribution 
# function (cdf) of $jobNumber to a cdf of a chi-square with 4 d.f. (as expected 
# for a distribution with four Gaussian variables.
#
# See the JKS paper:
# P. Jaranowski, A. Krolak, and B. F. Schutz, Phys. Rev. D 58, 063001 (1998).
#
#------------------------------------------------------------------------------

from __future__ import print_function
from sys import argv
from scipy import stats
from scipy.linalg import toeplitz
from scipy.special import gammaln
from math import log,pow,sqrt,exp,ceil,gamma
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import csv

# ---------- 1: Read in job  ----------------------------------------------------

jobNumber = str(argv[1])
jobFileName = "search_histogram.txt." + jobNumber  # Read in jobNumber from command line

# Note the following uses os.walk() in a way which will search for the first match 
# from the current working directory. This may cause unthrown errors!  

jobFullName = []
for walkRoot, walkDirs, walkFiles in os.walk(os.getcwd() + "/jobs/search/" + str(int(np.ceil(int(jobNumber)/250))) ):
    if jobFileName in walkFiles:
        jobFullName = os.path.join(walkRoot, jobFileName)


if not jobFullName:
    print("Error: could not find file ", jobFileName)
    exit(1)


# ---------- 2: Load in search histogram for job --------------------------------

# Get column data:									    

xHist = []
fHist = []

##############################################################
# This would be the proper way to index the x-values, i.e.
# as a fixed-point number to 3 d.p.
# This would use the 'decimal' package but this isn't defined well
# in SciPy and is thus incompatible with the scipy.stats package.
# I've left this here in comments just in case.
# R.I. 20140111
#
#from decimal import *
#getcontext().prec = 3  # Set the precision for fixed point
#with open(jobFullName, 'rb') as input:	
#	for row in csv.reader(input, delimiter=' '): 
#		if not row[0][0] == '%':
#			xHist.append(Decimal(row[0]))
#			fHist.append(int(row[2]))
#
#################################################################

with open(jobFullName, 'rb') as input:
    for row in csv.reader(input, delimiter=' '):
        if not row[0][0] == '%':
            xHist.append(float(row[0]))
            fHist.append(int(row[2]))


# ---------- 3: Calculate cdf of job  -------------------------------------------

dataDist = []
dataDist.append(fHist[0])
for xLabel in range(1, len(xHist)):
	#dataDist.append(dataDist[xLabel - 1] + fHist[xLabel])
	dataDist.append(dataDist[xLabel - 1] + fHist[xLabel])

#### Normalise the data

maxDataDist = float(max(dataDist))  # Need float so division is floating point(!)
dataDistNrml = []
for a in dataDist:
	dataDistNrml.append(a/maxDataDist)


# ---------- 4: Calculate cdf of 4 d.f. Chi^2 for same data points --------------


nullDist = stats.chi2.cdf(xHist,4,loc=0,scale=1)

# ---------- 5: Take difference, find maximum -----------------------------------


ksPlot = abs(dataDistNrml - nullDist)   # This may be used to show where the KS stat is maximum 

ksStat = max(ksPlot)

# ---------- 6: Get K-S statistic -----------------------------------------------

# TODO Check that ksStat actually is the K-S stat!!!

Ntot = len(xHist)
def criticalValue(alpha, Ntot):
    """Get the critical value from the K-S stat.
    This is valid for Ntot > 20 (Miller's approximation).
    This comes from the Matlab implementation of
    the (one-sided) K-S test.
    References:
    %   Massey, F.J., (1951) "The Kolmogorov-Smirnov Test for Goodness of Fit",
    %         Journal of the American Statistical Association, 46(253):68-78.
    %   Miller, L.H., (1956) "Table of Percentage Points of Kolmogorov Statistics",
    %         Journal of the American Statistical Association, 51(273):111-121.
    %   Marsaglia, G., W.W. Tsang, and J. Wang (2003), "Evaluating Kolmogorov`s
    %         Distribution", Journal of Statistical Software, vol. 8, issue 18.
    Note that for two-sided tests, you need to halve alpha.  """
    A =  0.09037*pow(-1.0*log(alpha, 10), 1.5 ) + 0.01515*pow( log(alpha, 10), 2 ) - 0.08467*alpha - 0.11143
    asymptoticStat =  sqrt(-0.5*log(alpha)/Ntot)
    criticalValue  =  asymptoticStat - 0.16693/Ntot - pow(A/Ntot,1.5)
    return min(criticalValue , 1 - alpha)

def noiseProb(ksStat, Ntot):
    """pvalue = noiseProb(ksStat, Ntot)
    Returns probability of being in same distribution for given K-S statistic"""
    s = Ntot*pow(ksStat, 2)
    # For d values that are in the far tail of the distribution (i.e.
    # p-values > .999), the following lines will speed up the computation
    # significantly, and provide accuracy up to 7 digits.
    if (s > 7.24) or ( (s > 3.76) and (Ntot > 99) ):
        return 2.0*exp( -1.0 * (2.000071 + 0.331/sqrt(Ntot) + 1.409/Ntot )*s )
    else:
        # Express d as d = (k-h)/Ntot, where k is a +ve integer and 0 < h < 1.
        k = ceil(ksStat*Ntot)
        h = k - ksStat*Ntot
        m = 2*k - 1.0
        # Create the H matrix, which describes the CDF, as described in Marsaglia,
        # et al. 
        if m > 1.0:
            c = [1.0/gamma(x) for x in np.linspace(2,m+1,m)]
            r = np.linspace(0,0,m)
            r[0] = 1.0
            r[1] = 1.0
            T = toeplitz(c,r)
            T[:,0] = T[:,0] - [ pow(h, x)/gamma(x + 1) for x in np.linspace(1,m,m) ]
            #T[m-1, :] = np.fliplr( T[:,0] )
            # Try this because fliplr is throwing a fit...
            L = len(T[:,0])
            T[m-1, :] = [ T[L-x-1,0] for x in range(L) ]
            T[m-1, 1] = (1.0 - 2.0*pow(h,m) + pow( max(0, 2.0*h - 1.0), m ) )/gamma(m+1)
        else:
            T = (1.0 - 2.0*pow( h, m ) + pow( max(0,2.0*h - 1), m) )/gamma(m+1)
            T=np.array([T,T])  # Make T appear to be a 2D array
        # Scaling before raising the matrix to a power
        if not np.isscalar(T):
            lmax = abs( np.amax(np.linalg.eig(T)[0]) )
            T = np.power( T/lmax, Ntot )
        else:
            lmax = 1.0
        # Pr(DNtot < d) = Ntot!/Ntot * tkk ,  where tkk is the kth element of TNtot = T^Ntot.
        # p-value = Pr(Dn > d) = 1-Pr(Dn < d)
        return 1.0 - exp(gammaln(Ntot+1) + Ntot*log(lmax) - Ntot*log(Ntot)) * T[k-1,k-1]



##---------------Test mode only ------------
#print( "The K-S statistic for job {%s is %f}".format(str(argv[1]), ksStat )

print("The K-S statistic for job " + jobNumber + " is: " + str(ksStat))


##----------------------------------------------


#############################################################################################
#  Plot/save results if prompted   #########
############################################
# Note that, without using the matplotlib API in a more detailed way,
# The following requires the $DISPLAY variable to be set. This is
# *not* a default setting when logging into remote machines!!!
# You can set this using e.g. X11 tunnelling (i.e. ssh -Y user@hostname )
###################################################################################


# TODO in light of having to set $DISPLAY , plotting should be optional.


plt.figure(1)
plt.subplot(211)   # Delete the subplot references to keep single plot mode

plt.plot(xHist, nullDist, "r-", label="Null distribution") 
plt.plot(xHist, dataDistNrml, "k--", label="Data")

plt.axis([-(max(xHist)-min(xHist))*0.05, max(xHist)*1.05, -0.01, max(dataDistNrml)*1.10])


xForPlot = np.linspace(min(xHist), max(xHist), 5)  # Make 5 marks on abscissa and ordinate 
yForPlot = np.linspace(0, max(dataDistNrml), 5)
x2DecPlcs = ['%.2f' % a for a in xForPlot ]
y2DecPlcs = ['%.2f' % a for a in yForPlot ]
plt.xticks(xForPlot, x2DecPlcs) 
plt.yticks(yForPlot, y2DecPlcs)

plt.title("Calculated cdf for $\chi ^2 _4 $ distribution of job {0}".format(jobNumber))
#plt.xlabel("2F value")
plt.ylabel("Cumulative frequency")


## Add a legend 
legend = plt.legend(loc='center right', shadow=True)

#frame = legend.get_frame()   # Some probably overly sophisticated additions to the legend
#frame.set_facecolor('0.90')
## Set the fontsize
#for label in legend.get_texts():
#	    label.set_fontsize('large')
#for label in legend.get_lines():
#	label.set_linewidth(1.5)  # the legend line width

#plt.draw()
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.grid(True)
#plt.show()

plt.text( 0.6*max(xHist), 0.2*max(dataDistNrml), 'K-S statistic: %.4f'%ksStat)    # Add the KS statistic on the plot

###########################################
### The following is for a second subplot showing where/how the KS stat gets large
### To remove, comment out, and remove subplot references above
plt.subplot(212)
plt.plot(xHist, ksPlot, 'r-')

plt.axis([-(max(xHist)-min(xHist))*0.05, max(xHist)*1.05, -(max(ksPlot)-min(ksPlot))*0.05, max(ksPlot)*1.10])

yForPlotKS = np.linspace(min(ksPlot), max(ksPlot), 5)
y3DecPlcsKS = ['%.3f' % a for a in yForPlotKS ]
plt.xticks(xForPlot, x2DecPlcs)
plt.yticks(yForPlotKS, y3DecPlcsKS)

#plt.title("K-S statistic"))
plt.xlabel("2F value")
plt.ylabel("|Difference|")

###########################################################################################



plt.draw()
plt.savefig("ksStat_" + jobNumber + ".png", dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format="png", transparent=False, bbox_inches=None, pad_inches=0.5, frameon=None)

#
#############################################################################################

exit(0)

# End of ksTest.py
