# fits1

This exercise will begin to introduce non linear data fitting.  The examples will be based on the ROOT libraries.  Similar functionality for everything we see can be found in the numpy + scipy + lmfit + matplotlib modules. But the examples we will follow over the next few projects offer several advantages:
* far more natural histogramming tools
* completely automated fitting options ("one liners")
* convenient methods for defining custom functions and controlling the convergence of the fit
* detailed and consistent documentation
* and a low level interface to modify the objective function, running in fully optimized compiled code

You are welcome to modify the provided code for your projects and to use other packages.  Where applicable alternate examples will be included. 

* **fit1.C**: C++ function to generate random data according to a normal distribution with mean=20, sigma=10. <br> A fit is performed using the built in Gaussian model in ROOT.  Then the parameter values, their uncertainteis, and the p-value for the fit are extracted.  To run this code type ```root fit1.C``` or if you are already running ROOT, type ```.X fit1.C```  
* **fit1.py**: The same code using the python interface, run this example using ```python fit1.C```.
* For a contrast see **fit1mpl.py** for a version using matplotlib+scipy.  
* readhist.C(py):  Examples for reading the histogram files given in this example 
* ParamUnceratinties.ipynb : a guided tutorial towards most of what you will be coding in this week's exercise.
* LLexample.ipynb : a notebook giving an example for calculating (N)LLs
* TH1hist2Numpy.ipynb : an example for converting a ROOT histogram to numpy arrays

Note that from ROOT you can type ```new TBrowser()``` or in Python r.TBrowser() to get a graphical browser that allows you to look at what's contained in the TFiles.



STUDENT ANSWERS TO CANVAS QUESTIONS:

EXERCISE 1:
Q: How does the 1 sigma width of this distribution compare to the typical size of the uncertainty reported for this fit parameter?
A: Comparing our plot statistics with our printed uncertainties for the mean, they are the same. This seems reasonable, because I would expect that if the mean value has an uncertainty of - in this case - 0.369, then most mean values should only differ by that much, which is nothing else but one sigma in the distribution of mean values.

EXERCISE 2:
Q: How do your results compare to the expected values in each case? How do the distributions of the parameter values from the fits compare to the estimated uncertainty on the fit parameters?
A: As anticipated, the NLL method remained more consistently close to the expected value with a relatively small uncertainty, whereas the chi-squared method still centered on the correct value, but had far more variation and much higher errors from experiment to experiment.

EXERCISE 3:
Q: Estimate the "p-value" of your original fit.
A: Run to run, our p-value hovers around 0.5

Q: Give a comparison of the errors reported by the fitter and the results of your scans.
A: Both NLL and Chi2 did remarkably well at recapturing the ROOT fit mean values from minimizing the contour. Specifically, we found for NLL, the ROOT fit gave mu = 53.420 +/- 2.088 and the contour was at mu = 53.315. The chi-squared ROOT fit gave mu = 49.781 +/- 0.323 and the contour gave mu = 49.797.



