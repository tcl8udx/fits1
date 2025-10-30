import ROOT as r
import numpy as np

def fit1(entries=10, ntrials=1000, save=True):
    # entries is the number of random samples filled into the histogram

    randomHist1 = r.TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100)
    randomHist1.Sumw2()
    generator=r.TRandom2(0)  # parameter == seed, 0->use clock

    # Arrays to store results from each trial
    #chi2_values = []
    #prob_values = []
    #reduced_chi2_values = []
    chi2_mean_values = []
    chi2_mean_errors = []
    NLL_mean_values = []
    NLL_mean_errors = []
    
    print(f"Running {ntrials} experiments with {entries} points each...")
    print("=="*40)
    print(f"{'Trial':<10}{'Chi^2 Mean':<15}{'Likelihood Mean':<15}")

    for j in range(ntrials):
        randomHist1.Reset();   # reset histogram bin content to 0
        for i in range(entries):
            randomHist1.Fill(generator.Gaus(50,10)) # params: mean, sigma
            
        randomHist1.Fit("gaus", "Q") # Chi^2 default
        fitfunc_chi2 = randomHist1.GetFunction("gaus")
        chi2_mean = fitfunc_chi2.GetParameter(1)  # mean is parameter 1
        chi2_mean_error = fitfunc_chi2.GetParError(1)

        randomHist1.Fit("gaus", "QL") # NLL fit
        fitfunc_NLL = randomHist1.GetFunction("gaus")
        NLL_mean = fitfunc_NLL.GetParameter(1)  # mean is parameter 1
        NLL_mean_error = fitfunc_NLL.GetParError(1)
        
    
        if fitfunc_chi2:  # Check if fit was successful
            # Store results
            chi2_mean_values.append(chi2_mean)
            chi2_mean_errors.append(chi2_mean_error)

        if fitfunc_NLL:  # Check if fit was successful
            # Store results
            NLL_mean_values.append(NLL_mean)
            NLL_mean_errors.append(NLL_mean_error)
            
            # Print every 100th trial
            if j % 100 == 0:
                print(f"{j:<10}{chi2_mean:<15.2f}{NLL_mean:<15.3f}")

    print(f"\nCompleted {len(chi2_mean_values)} successful fits")
    
    if save:
        canvas = r.TCanvas("canvas", "Fit Results", 1200, 1200)
        canvas.Divide(2, 1)
        
        # Subplot 1: Fitted chi2 mean distribution
        canvas.cd(1)
        chi2_h_mean = r.TH1F("chi2_h_mean", "Fitted Chi2 Mean Distribution;Fitted Mean;Frequency", 50, 20, 80)
        for val in chi2_mean_values:
            chi2_h_mean.Fill(val)
        chi2_h_mean.Draw()
        chi2_h_mean.Fit("gaus", "Q")
        
        # Add statistics
        #stats1 = r.TPaveText(0.6, 0.6, 0.89, 0.89, "NDC")
        #stats1.AddText(f"Mean: {np.mean(reduced_chi2_values):.3f}")
        #stats1.AddText(f"Std: {np.std(reduced_chi2_values):.3f}")
        #stats1.AddText(f"Entries: {len(reduced_chi2_values)}")
        #stats1.SetFillColor(0)
        #stats1.Draw()
        
        # Subplot 2: Fitted NLL mean distribution
        canvas.cd(2)
        NLL_h_mean = r.TH1F("NLL_h_mean", "Fitted NLL Mean Distribution;Fitted Mean;Frequency", 50, 30, 70)
        for val in NLL_mean_values:
            NLL_h_mean.Fill(val)
        NLL_h_mean.Draw()
        NLL_h_mean.Fit("gaus", "Q")
        
        #stats2 = r.TPaveText(0.6, 0.6, 0.89, 0.89, "NDC")
        #stats2.AddText(f"Mean: {np.mean(prob_values):.3f}")
        #stats2.AddText(f"Std: {np.std(prob_values):.3f}")
        #stats2.AddText(f"Entries: {len(prob_values)}")
        #stats2.SetFillColor(0)
        #stats2.Draw()
        
        
        # Save to PDF
        canvas.SaveAs("result2.pdf")
        print(f"\nResults saved to result1.pdf")
        
        # Print summary statistics
        print("\n" + "=="*40)
        print("SUMMARY STATISTICS:")
        print("=="*40)
        print(f"Fitted Chi2 Mean: {np.mean(chi2_mean_values):.3f} +/- {np.std(chi2_mean_values):.3f}")
        print(f"Fitted Likelihood Mean: {np.mean(NLL_mean_values):.3f} +/- {np.std(NLL_mean_values):.3f}")
        print(f"Expected Mean: 50.0")
        print(f"Chi2 Mean Error: {np.mean(chi2_mean_errors):.3f} +/- {np.std(chi2_mean_errors):.3f}")
        print(f"NLL Mean Error: {np.mean(NLL_mean_errors):.3f} +/- {np.std(NLL_mean_errors):.3f}")
        print(f"Std Error (theory): {10.0/np.sqrt(entries):.3f}")
        
        return canvas, (chi2_h_mean, NLL_h_mean)
    
    return None
    
# **************************************

if __name__ == "__main__":
    result = fit1()
    input("hit Enter to exit")


# example of plotting a similar histogram with error bars using numpy/matplotlib
# then using lmfit to perform the fit

#import numpy as np
#from matplotlib import pyplot as plt
#
#entries=1000
#
#vals=np.random.normal(loc=50, scale=10, size=entries)
#y,binEdges=np.histogram(vals, bins=50, range=(1,100))
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#ey         = np.sqrt(y)
#width      = binEdges[1]-binEdges[0]
#plt.bar(bincenters, y, width=width, color='r', yerr=ey)
#plt.show(block=False)
