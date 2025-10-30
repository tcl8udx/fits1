import ROOT as r
import numpy as np

def fit1(entries=1000, ntrials=1000, save=True):
    # entries is the number of random samples filled into the histogram

    randomHist1 = r.TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100)
    randomHist1.Sumw2()
    generator=r.TRandom2(0)  # parameter == seed, 0->use clock

    # Arrays to store results from each trial
    chi2_values = []
    prob_values = []
    reduced_chi2_values = []
    mean_values = []
    mean_errors = []
    
    print(f"Running {ntrials} experiments with {entries} points each...")
    print("=="*40)
    print(f"{'Trial':<10}{'Chi^2':<15}{'Probability':<15}{'Reduced Chi^2':<15}{'Mean':<15}")

    for j in range(ntrials):
        randomHist1.Reset();   # reset histogram bin content to 0
        for i in range(entries):
            randomHist1.Fill(generator.Gaus(50,10)) # params: mean, sigma
        randomHist1.Fit("gaus","Q") # Q = quiet mode

        fitfunc = randomHist1.GetFunction("gaus")
    
        if fitfunc:  # Check if fit was successful
            chi2 = fitfunc.GetChisquare()
            ndf = fitfunc.GetNDF()
            prob = fitfunc.GetProb()
            reduced_chi2 = chi2 / ndf if ndf > 0 else 0
            mean = fitfunc.GetParameter(1)  # mean is parameter 1
            mean_error = fitfunc.GetParError(1)
            
            # Store results
            chi2_values.append(chi2)
            prob_values.append(prob)
            reduced_chi2_values.append(reduced_chi2)
            mean_values.append(mean)
            mean_errors.append(mean_error)
            
            # Print every 100th trial
            if j % 100 == 0:
                print(f"{j:<10}{chi2:<15.2f}{prob:<15.3f}{reduced_chi2:<15.2f}{mean:<15.2f}")

    print(f"\nCompleted {len(reduced_chi2_values)} successful fits")
    
    if save and len(reduced_chi2_values) > 0:
        canvas = r.TCanvas("canvas", "Fit Results", 1200, 1200)
        canvas.Divide(2, 2)
        
        # Subplot 1: Reduced Chi-squared distribution
        canvas.cd(1)
        h_redchi2 = r.TH1F("h_redchi2", "Reduced #chi^{2} Distribution;Reduced #chi^{2};Frequency", 50, 0, 3)
        for val in reduced_chi2_values:
            h_redchi2.Fill(val)
        h_redchi2.Draw()
        
        # Add statistics
        #stats1 = r.TPaveText(0.6, 0.6, 0.89, 0.89, "NDC")
        #stats1.AddText(f"Mean: {np.mean(reduced_chi2_values):.3f}")
        #stats1.AddText(f"Std: {np.std(reduced_chi2_values):.3f}")
        #stats1.AddText(f"Entries: {len(reduced_chi2_values)}")
        #stats1.SetFillColor(0)
        #stats1.Draw()
        
        # Subplot 2: Chi-squared probability distribution
        canvas.cd(3)
        h_prob = r.TH1F("h_prob", "#chi^{2} Probability Distribution;Probability;Frequency", 50, 0, 1)
        for val in prob_values:
            h_prob.Fill(val)
        h_prob.Draw()
        
        #stats2 = r.TPaveText(0.6, 0.6, 0.89, 0.89, "NDC")
        #stats2.AddText(f"Mean: {np.mean(prob_values):.3f}")
        #stats2.AddText(f"Std: {np.std(prob_values):.3f}")
        #stats2.AddText(f"Entries: {len(prob_values)}")
        #stats2.SetFillColor(0)
        #stats2.Draw()
        
        # Subplot 3: Fitted mean distribution
        canvas.cd(2)
        h_mean = r.TH1F("h_mean", "Fitted Mean Distribution;Fitted Mean;Frequency", 50, 48, 52)
        for val in mean_values:
            h_mean.Fill(val)
        h_mean.Draw()
        h_mean.Fit("gaus", "Q")
        
        #stats3 = r.TPaveText(0.6, 0.6, 0.89, 0.89, "NDC")
        #stats3.AddText(f"Mean: {np.mean(mean_values):.3f}")
        #stats3.AddText(f"Std: {np.std(mean_values):.3f}")
        #stats3.AddText(f"True Mean: 50.0")
        #stats3.SetFillColor(0)
        #stats3.Draw()
        
        # Subplot 4: Error on mean distribution
        canvas.cd(4)
        h_error = r.TH1F("h_error", "Error on Mean Distribution;Error on Mean;Frequency", 50, 0, 1)
        for val in mean_errors:
            h_error.Fill(val)
        h_error.Draw()
        
        #stats4 = r.TPaveText(0.6, 0.6, 0.89, 0.89, "NDC")
        #stats4.AddText(f"Mean: {np.mean(mean_errors):.3f}")
        #stats4.AddText(f"Std: {np.std(mean_errors):.3f}")
        #stats4.AddText(f"Entries: {len(mean_errors)}")
        #stats4.SetFillColor(0)
        #stats4.Draw()
        
        # Save to PDF
        canvas.SaveAs("result1.pdf")
        print(f"\nResults saved to result1.pdf")
        
        # Print summary statistics
        print("\n" + "=="*40)
        print("SUMMARY STATISTICS:")
        print("=="*40)
        print(f"Reduced Chi^2: {np.mean(reduced_chi2_values):.3f} +/- {np.std(reduced_chi2_values):.3f}")
        print(f"Chi^2 Probability: {np.mean(prob_values):.3f} +/- {np.std(prob_values):.3f}")
        print(f"Fitted Mean: {np.mean(mean_values):.3f} +/- {np.std(mean_values):.3f}")
        print(f"Expected Mean: 50.0")
        print(f"Mean Error: {np.mean(mean_errors):.3f} +/- {np.std(mean_errors):.3f}")
        print(f"Std Dev of Means: {np.std(mean_values):.3f}")
        print(f"Std Error (theory): {10.0/np.sqrt(entries):.3f}")
        
        return canvas, (h_redchi2, h_prob, h_mean, h_error)
    
    return None
    
# **************************************

if __name__ == "__main__":
    result = fit1(entries=1000, ntrials=1000, save=True)
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
