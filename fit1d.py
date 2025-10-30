import ROOT as r
import numpy as np

def calculate_nll(histogram, mean, sigma, amplitude):
    """
    Calculate the Negative Log-Likelihood for a histogram given Gaussian parameters.
    Returns -2*NLL.
    """
    nll = 0.0
    nbins = histogram.GetNbinsX()
    model = r.TF1("model", "gaus", 0,)
    model.SetParameter(0, amplitude)
    model.SetParameter(1, mean)
    model.SetParameter(2, sigma)
    
    for i in range(1, nbins + 1):
        bin_content = histogram.GetBinContent(i)
        bin_center = histogram.GetBinCenter(i)
        nu = model.Eval(bin_center)
        
        if nu <= 0:
            continue
            
        if bin_content > 0:
            # Poisson likelihood: -ln(L) = sum[expected - observed*ln(expected)]
            nll += 2*(nu - bin_content*np.log(nu))

        else:
            nll += 2*nu
    
    return nll

def calculate_chi2(histogram, mean, sigma, amplitude):
    """
    Calculate Chi-Squared for a histogram given Gaussian parameters.
    """
    chi2 = 0.0
    nbins = histogram.GetNbinsX()
    
    for i in range(1, nbins + 1):
        bin_content = histogram.GetBinContent(i)
        bin_error = histogram.GetBinError(i)
        bin_center = histogram.GetBinCenter(i)
        
        if bin_error > 0:  # Only include bins with valid errors
            # Calculate expected value from Gaussian
            expected = amplitude * np.exp(-0.5 * ((bin_center - mean) / sigma)**2)
            
            # Chi-squared = sum[(observed - expected)^2 / error^2]
            chi2 += ((bin_content - expected) / bin_error)**2
    
    return chi2

def plot_likelihood_contour(rootfile_path, histogram_name, save=True):
    """
    Plot -2*ln(L) contour around minimum by varying the mean parameter.
    """
    print(f"="*60)
    print("LIKELIHOOD CONTOUR ANALYSIS")
    print(f"="*60)
    
    # Load histogram
    print(f"\nLoading histogram from {rootfile_path}...")
    file = r.TFile(rootfile_path)
    hist = file.Get(histogram_name)
    
    if not hist:
        print(f"Error: Could not find histogram '{histogram_name}'")
        return None
    
    hist = hist.Clone("hist_clone")
    hist.SetDirectory(0)
    hist.Sumw2()
    
    # Perform NLL fit
    print("Performing NLL fit...")
    hist.Fit("gaus", "QL")
    fitfunc = hist.GetFunction("gaus")
    
    if not fitfunc:
        print("Error: Fit failed")
        return None
    
    best_amp = fitfunc.GetParameter(0)
    best_mean = fitfunc.GetParameter(1)
    best_sigma = fitfunc.GetParameter(2)
    mean_error = fitfunc.GetParError(1)
    
    print(f"\nBest fit parameters:")
    print(f"  Amplitude: {best_amp:.2f}")
    print(f"  Mean: {best_mean:.3f} +/- {mean_error:.3f}")
    print(f"  Sigma: {best_sigma:.3f}")
    
    # Calculate NLL at minimum
    nll_min = calculate_nll(hist, best_mean, best_sigma, best_amp)
    minus_2lnL_min = nll_min
    
    print(f"-2*ln(L) at minimum: {minus_2lnL_min:.2f}")
    
    # Scan mean parameter around best fit value
    # Scan range: typically +/- 3-5 sigma to capture up to -2lnL_min + 4
    scan_range = 5 * mean_error
    n_points = 100
    mean_values = np.linspace(best_mean - scan_range, best_mean + scan_range, n_points)
    minus_2lnL_values = []
    
    print(f"\nScanning mean from {mean_values[0]:.3f} to {mean_values[-1]:.3f}...")
    
    for mean_val in mean_values:
        nll = calculate_nll(hist, mean_val, best_sigma, best_amp)
        minus_2lnL_values.append(nll)

    # Find mean at minimum
    imin = np.argmin(minus_2lnL_values)
    mu_best_scan = mean_values[imin]
    print(f"Mean value at NLL contour minimum: {mu_best_scan}")
    
    # Create plot
    canvas = r.TCanvas("canvas_likelihood", "Likelihood Contour", 800, 600)
    
    # Create TGraph
    graph = r.TGraph(n_points, mean_values, np.array(minus_2lnL_values))
    graph.SetTitle("-2ln(L) vs Mean Parameter;Mean;-2ln(L)")
    graph.SetLineColor(r.kBlue)
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(0)
    
    # Find y-axis range
    y_min = minus_2lnL_min - 0.5
    y_max = minus_2lnL_min + 5
    graph.GetYaxis().SetRangeUser(y_min, y_max)
    
    graph.Draw("AL")
    
    # Draw horizontal lines at -2lnL_min + 1 and -2lnL_min + 4
    x_min = mean_values[0]
    x_max = mean_values[-1]
    
    line_min = r.TLine(x_min, minus_2lnL_min, x_max, minus_2lnL_min)
    line_min.SetLineColor(r.kRed)
    line_min.SetLineStyle(2)
    line_min.SetLineWidth(2)
    line_min.Draw("same")
    
    line_1sigma = r.TLine(x_min, minus_2lnL_min + 1, x_max, minus_2lnL_min + 1)
    line_1sigma.SetLineColor(r.kGreen + 2)
    line_1sigma.SetLineStyle(2)
    line_1sigma.SetLineWidth(2)
    line_1sigma.Draw("same")
    
    line_2sigma = r.TLine(x_min, minus_2lnL_min + 4, x_max, minus_2lnL_min + 4)
    line_2sigma.SetLineColor(r.kOrange + 1)
    line_2sigma.SetLineStyle(2)
    line_2sigma.SetLineWidth(2)
    line_2sigma.Draw("same")
    
    # Add legend
    legend = r.TLegend(0.55, 0.65, 0.89, 0.89)
    legend.AddEntry(graph, "-2ln(L)", "l")
    legend.AddEntry(line_min, "Minimum", "l")
    legend.AddEntry(line_1sigma, "#Delta(-2ln(L)) = 1 (68% CL)", "l")
    legend.AddEntry(line_2sigma, "#Delta(-2ln(L)) = 4 (95% CL)", "l")
    legend.Draw()
    
    # Add text box with best fit value
    text = r.TPaveText(0.15, 0.75, 0.45, 0.89, "NDC")
    text.AddText(f"Best fit mean: {best_mean:.3f}")
    text.AddText(f"-2ln(L)_{{min}}: {minus_2lnL_min:.2f}")
    text.SetFillColor(0)
    text.SetTextAlign(12)
    text.Draw()
    
    if save:
        canvas.SaveAs("likelihood_contour.pdf")
        print(f"\nLikelihood contour saved to likelihood_contour.pdf")
    
    return canvas, graph, hist

def plot_chi2_contour(rootfile_path, histogram_name, save=True):
    """
    Plot Chi-Squared contour around minimum by varying the mean parameter.
    """
    print(f"\n{'='*60}")
    print("CHI-SQUARED CONTOUR ANALYSIS")
    print(f"{'='*60}")
    
    # Load histogram
    print(f"\nLoading histogram from {rootfile_path}...")
    file = r.TFile(rootfile_path)
    hist = file.Get(histogram_name)
    
    if not hist:
        print(f"Error: Could not find histogram '{histogram_name}'")
        return None
    
    hist = hist.Clone("hist_clone_chi2")
    hist.SetDirectory(0)
    hist.Sumw2()
    
    # Perform Chi-Squared fit
    print("Performing Chi-Squared fit...")
    hist.Fit("gaus", "Q")
    fitfunc = hist.GetFunction("gaus")
    
    if not fitfunc:
        print("Error: Fit failed")
        return None
    
    best_amp = fitfunc.GetParameter(0)
    best_mean = fitfunc.GetParameter(1)
    best_sigma = fitfunc.GetParameter(2)
    mean_error = fitfunc.GetParError(1)
    
    print(f"\nBest fit parameters:")
    print(f"  Amplitude: {best_amp:.2f}")
    print(f"  Mean: {best_mean:.3f} +/- {mean_error:.3f}")
    print(f"  Sigma: {best_sigma:.3f}")
    
    # Calculate Chi2 at minimum
    chi2_min = calculate_chi2(hist, best_mean, best_sigma, best_amp)
    
    print(f"\nChi-Squared at minimum: {chi2_min:.2f}")
    
    # Scan mean parameter
    scan_range = 5 * mean_error
    n_points = 100
    mean_values = np.linspace(best_mean - scan_range, best_mean + scan_range, n_points)
    chi2_values = []
    
    print(f"\nScanning mean from {mean_values[0]:.3f} to {mean_values[-1]:.3f}...")
    
    for mean_val in mean_values:
        chi2 = calculate_chi2(hist, mean_val, best_sigma, best_amp)
        chi2_values.append(chi2)

    # Find mean at minimum
    imin = np.argmin(chi2_values)
    mu_best_scan = mean_values[imin]
    print(f"Mean value at chi2 contour minimum: {mu_best_scan}")
    
    # Create plot
    canvas = r.TCanvas("canvas_chi2", "Chi-Squared Contour", 800, 600)
    
    # Create TGraph
    graph = r.TGraph(n_points, mean_values, np.array(chi2_values))
    graph.SetTitle("#chi^{2} vs Mean Parameter;Mean;#chi^{2}")
    graph.SetLineColor(r.kBlue)
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(0)
    
    # Find y-axis range
    y_min = chi2_min - 0.5
    y_max = chi2_min + 5
    graph.GetYaxis().SetRangeUser(y_min, y_max)
    
    graph.Draw("AL")
    
    # Draw horizontal lines
    x_min = mean_values[0]
    x_max = mean_values[-1]
    
    line_min = r.TLine(x_min, chi2_min, x_max, chi2_min)
    line_min.SetLineColor(r.kRed)
    line_min.SetLineStyle(2)
    line_min.SetLineWidth(2)
    line_min.Draw("same")
    
    line_1sigma = r.TLine(x_min, chi2_min + 1, x_max, chi2_min + 1)
    line_1sigma.SetLineColor(r.kGreen + 2)
    line_1sigma.SetLineStyle(2)
    line_1sigma.SetLineWidth(2)
    line_1sigma.Draw("same")
    
    line_2sigma = r.TLine(x_min, chi2_min + 4, x_max, chi2_min + 4)
    line_2sigma.SetLineColor(r.kOrange + 1)
    line_2sigma.SetLineStyle(2)
    line_2sigma.SetLineWidth(2)
    line_2sigma.Draw("same")
    
    # Add legend
    legend = r.TLegend(0.55, 0.65, 0.89, 0.89)
    legend.AddEntry(graph, "#chi^{2}", "l")
    legend.AddEntry(line_min, "Minimum", "l")
    legend.AddEntry(line_1sigma, "#Delta#chi^{2} = 1 (68% CL)", "l")
    legend.AddEntry(line_2sigma, "#Delta#chi^{2} = 4 (95% CL)", "l")
    legend.Draw()
    
    # Add text box
    text = r.TPaveText(0.15, 0.75, 0.45, 0.89, "NDC")
    text.AddText(f"Best fit mean: {best_mean:.3f}")
    text.AddText(f"#chi^{{2}}_{{min}}: {chi2_min:.2f}")
    text.SetFillColor(0)
    text.SetTextAlign(12)
    text.Draw()
    
    if save:
        canvas.SaveAs("chi2_contour.pdf")
        print(f"\nChi-squared contour saved to chi2_contour.pdf")
    
    return canvas, graph, hist

# **************************************
if __name__ == "__main__":
    # Plot likelihood contour for histo25.root
    print("Processing likelihood contour...")
    result_likelihood = plot_likelihood_contour("histo25.root", "randomHist1", save=False)
    
    print("\n")
    
    # Plot chi-squared contour for histo1k.root
    print("Processing chi-squared contour...")
    result_chi2 = plot_chi2_contour("histo1k.root", "randomHist1", save=False)

    # Save both to a single multi-page PDF
    result_likelihood[0].Print("result4.pdf(")  # Open PDF with first page
    result_chi2[0].Print("result4.pdf)")        # Close PDF with second page
    
    input("\nHit Enter to exit")