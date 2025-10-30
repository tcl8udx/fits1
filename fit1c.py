import ROOT as r
import numpy as np

def calculate_nll(histogram, mean, sigma, amplitude):
    """
    Calculate the Negative Log-Likelihood for a histogram given Gaussian parameters.
    
    NLL = -sum_i[ln(P_i)] where P_i is the probability for bin i
    For a Gaussian: P_i = (N/sqrt(2*pi*sigma^2)) * exp(-(x_i - mu)^2 / (2*sigma^2)) * bin_width
    
    Parameters:
    - histogram: ROOT TH1F histogram
    - mean: Gaussian mean parameter
    - sigma: Gaussian sigma parameter  
    - amplitude: Gaussian amplitude (normalization)
    """
    nll = 0.0
    nbins = histogram.GetNbinsX()
    bin_width = histogram.GetBinWidth(1)
    
    for i in range(1, nbins + 1):  # ROOT bins start at 1
        bin_content = histogram.GetBinContent(i)
        bin_center = histogram.GetBinCenter(i)
        
        if bin_content > 0:  # Only include bins with data
            # Calculate expected value from Gaussian
            expected = amplitude * np.exp(-0.5 * ((bin_center - mean) / sigma)**2)
            
            # Add small epsilon to avoid log(0)
            if expected > 0:
                # Poisson likelihood: -ln(L) = sum[expected - observed*ln(expected) + ln(observed!)]
                # Ignoring the ln(observed!) constant term:
                nll += expected - bin_content * np.log(expected)
    
    return nll


def nll_consistency_test(rootfile_path, histogram_name, ntrials=1000, save=True):
    """
    Test NLL fit consistency using pseudo-experiments.
    
    Parameters:
    - rootfile_path: path to .root file containing histogram
    - histogram_name: name of histogram in file
    - ntrials: number of pseudo-experiments
    - save: whether to save result plot
    """
    
    # Load the parent histogram
    print(f"Loading histogram from {rootfile_path}...")
    file = r.TFile(rootfile_path)
    parent_hist = file.Get(histogram_name)
    
    if not parent_hist:
        print(f"Error: Could not find histogram '{histogram_name}' in file")
        return None
    
    # Make a copy so we can close the file
    parent_hist = parent_hist.Clone("parent_hist_clone")
    parent_hist.SetDirectory(0)  # Detach from file
    
    # Get histogram properties
    nbins = parent_hist.GetNbinsX()
    xmin = parent_hist.GetXaxis().GetXmin()
    xmax = parent_hist.GetXaxis().GetXmax()
    total_entries = parent_hist.GetEntries()
    
    print(f"Histogram has {total_entries:.0f} entries in {nbins} bins from {xmin} to {xmax}")
    
    # Perform NLL fit on parent histogram
    print("\nPerforming NLL fit on parent histogram...")
    parent_hist.Fit("gaus", "QL")
    parent_fitfunc = parent_hist.GetFunction("gaus")
    
    if not parent_fitfunc:
        print("Error: Fit failed")
        return None
    
    parent_amp = parent_fitfunc.GetParameter(0)
    parent_mean = parent_fitfunc.GetParameter(1)
    parent_std = parent_fitfunc.GetParameter(2)
    
    print(f"Best fit parameters:")
    print(f"  Amplitude: {parent_amp:.2f}")
    print(f"  Mean: {parent_mean:.3f} +/- {parent_fitfunc.GetParError(1):.3f}")
    print(f"  Sigma: {parent_std:.3f} +/- {parent_fitfunc.GetParError(2):.3f}")
    
    # Calculate NLL for the parent histogram with best fit parameters
    parent_nll = calculate_nll(parent_hist, parent_mean, parent_std, parent_amp)
    print(f"\nNLL for parent histogram: {parent_nll:.2f}")
    
    # Generate pseudo-experiments
    print(f"\nGenerating {ntrials} pseudo-experiments...")
    generator = r.TRandom3(0)  # Use TRandom3 for better quality
    nll_values = []
    
    # Create histogram for pseudo-experiments
    pseudo_hist = r.TH1F("pseudo_hist", "Pseudo Histogram", nbins, xmin, xmax)
    
    for trial in range(ntrials):
        pseudo_hist.Reset()
        
        # Fill pseudo-histogram with random samples from best fit Gaussian
        for i in range(int(total_entries)):
            pseudo_hist.Fill(generator.Gaus(parent_mean, parent_std))
        
        # Calculate NLL for this pseudo-experiment using parent fit parameters
        pseudo_nll = calculate_nll(pseudo_hist, parent_mean, parent_std, parent_amp)
        nll_values.append(pseudo_nll)
        
        if trial % 100 == 0:
            print(f"  Trial {trial}: NLL = {pseudo_nll:.2f}")
    
    print(f"\nCompleted {len(nll_values)} pseudo-experiments")
    
    # Calculate statistics
    nll_array = np.array(nll_values)
    mean_nll = np.mean(nll_array)
    std_nll = np.std(nll_array)
    
    # Calculate p-value: fraction of pseudo-experiments with NLL >= parent NLL
    n_greater = np.sum(nll_array >= parent_nll)
    p_value = n_greater / len(nll_values)
    
    print("\n" + "=="*40)
    print("RESULTS:")
    print("=="*40)
    print(f"Parent NLL: {parent_nll:.2f}")
    print(f"Mean pseudo-experiment NLL: {mean_nll:.2f} +/- {std_nll:.2f}")
    print(f"P-value (fraction with NLL >= parent): {p_value:.3f}")
    print(f"Number of pseudo-experiments with NLL >= parent: {n_greater}/{len(nll_values)}")
    
    if p_value > 0.05:
        print(f"\nConclusion: Data is CONSISTENT with the fit (p = {p_value:.3f})")
    else:
        print(f"\nConclusion: Data may be INCONSISTENT with the fit (p = {p_value:.3f})")
    
    # Create plot
    if save:
        canvas = r.TCanvas("canvas", "NLL Distribution", 800, 600)
        
        # Create histogram of NLL values
        nll_min = min(nll_values)
        nll_max = max(nll_values)
        h_nll = r.TH1F("h_nll", "NLL Distribution from Pseudo-Experiments;NLL;Frequency", 
                       50, nll_min * 0.95, nll_max * 1.05)
        
        for val in nll_values:
            h_nll.Fill(val)
        
        h_nll.Draw()
        
        # Draw a line at the parent NLL value
        line = r.TLine(parent_nll, 0, parent_nll, h_nll.GetMaximum() * 1.1)
        line.SetLineColor(r.kRed)
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.Draw("same")
        
        # Add legend
        legend = r.TLegend(0.7, 0.5, 0.89, 0.69)
        legend.AddEntry(h_nll, "Pseudo-experiments", "f")
        legend.AddEntry(line, f"Parent NLL = {parent_nll:.1f}", "l")
        legend.SetFillColor(0)
        legend.Draw()
        
        # Add statistics box
        #stats = r.TPaveText(0.12, 0.7, 0.45, 0.89, "NDC")
        #stats.AddText(f"Mean NLL: {mean_nll:.2f}")
        #stats.AddText(f"Std NLL: {std_nll:.2f}")
        #stats.AddText(f"Entries: {len(nll_values)}")
        #stats.SetFillColor(0)
        #stats.SetTextAlign(12)
        #stats.Draw()
        
        canvas.SaveAs("result3.pdf")
        print(f"\nPlot saved to result3.pdf")
        
        return canvas, h_nll, parent_hist
    
    return None

# **************************************
if __name__ == "__main__":
    # Example usage
    result = nll_consistency_test("histo25.root", "randomHist1", ntrials=1000, save=True)
    input("Hit Enter to exit")