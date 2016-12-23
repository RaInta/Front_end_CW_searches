# Front_end_CW_searches
This is a front end to display the results of running the continuous gravitational wave searches.

The main output is a web page, using a really basic application of the Django web framework.

The root file is <a href="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/createHTML.py"><code>createHTML.py</code></a>, which calls a whole lot of other scripts. 

A couple of scripts, <a href="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/plotUpperLimits.py"><code>plotUpperLimits.py</code></a> and <a href="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/readOptimalStretch.py"><code>readOptimalStretch.py</code></a>  just read basic information to display as an overview. 

Another script, <a href="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/lookThresh.py"><code>lookThresh.py</code></a> is quite a sophisticated statistical exercise. 

It takes in the total number of templates, and computes (amongst other things) the maximum detection statistic ("twoF") expected for the search job. This would be relatively simple, as the twoF is a result of a maximum likelihood optimization over four (assumed i.i.d. random) variables, so a Chi^2 test with four degrees of freedom ought to give this result. However, these templates are not completely independent. So the effective dimensionality of independent templates is found by minimizing the Kolmogorov-Smirnoff (K-S) measure between a range of numbers of templates. If it's Gaussian noise, then the effective number of templates occurs at the minimum point of the K-S measure.    

To determine the Gaussianity of each search job, the Kolmogorov-Smirnoff (K-S) test for normailty is applied to each search job, using the <a href="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/ksTest.py"><code>ksTest.py</code></a> K-S test Python script.

An example of running <code>ksTest.py</code> on a job that is pure Gaussian noise can be seen here:

<img src="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/ksStat_P0_0.png" width="600">


For comparison, here is a search job that isn't Gaussian:

<img src="https://github.com/NotAFakeRa/Front_end_CW_searches/blob/master/ksStat_5374.png" width="600">


The significance of each of these jobs is worked out in the script (and output to stdout).

A real problem with the main script is that this was written while Django 1.7 was the latest thing. The way <a href="https://docs.djangoproject.com/en/1.8/topics/templates/">templates are handled was completely overhauled in 1.8 and later</a>, leading to backwards imcompatibility. However, all the other scripts are very useful.
