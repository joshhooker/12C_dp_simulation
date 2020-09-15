import array
import numpy as np
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error
import sys

from ROOT import gROOT, gStyle, TCanvas, TColor, TFile, TGaxis, TGraph, TGraphErrors, TH1, TLatex, TPad
from ROOT import kBlack, kBlue, kRed

def func(x, fit, sim):
    sim_copy = np.copy(sim)
    sim_copy /= x
    return mean_squared_error(fit, sim_copy)

def analysis(input_file, output_file):

    c1 = TCanvas('c1', 'Compare Simulation to Real Data', 200, 10, 700, 500)

    data_x, data_y, data_yerr = np.loadtxt("CQvalueTableShort.txt", skiprows=1, unpack=True)
    data_x = data_x.flatten('C'); data_x = data_x.astype('float')
    data_y = data_y.flatten('C'); data_y = data_y.astype('float')
    data_yerr = data_yerr.flatten('C'); data_yerr = data_yerr.astype('float')
    data_xerr = np.zeros(data_x.shape[0])
    data_xerr = data_xerr.flatten('C'); data_xerr = data_xerr.astype('float')

    try:
        input_tfile = TFile.Open(input_file, "read")
    except:
        return 100000.
    calcQValueHist = input_tfile.Get("calcQValue")
    i_size = calcQValueHist.GetSize()
    xaxis = calcQValueHist.GetXaxis()
    sim_x = []
    sim_y = []
    for i in range(i_size):
        binCenter = xaxis.GetBinCenter(i)
        binContent = calcQValueHist.GetBinContent(i)
        if (binCenter < -2.50):
            continue
        if (binCenter > 4.0):
            continue
        sim_x.append(binCenter)
        sim_y.append(binContent)

    sim_x = np.array(sim_x)
    sim_y = np.array(sim_y)

    print(data_x, len(data_x))
    print(sim_x, len(sim_x))

    result = minimize(func, 20, args=(data_y, sim_y), tol=1e-6)
    sim_y_scaled = np.copy(sim_y)/result.x[0]

    mse = mean_squared_error(data_y, sim_y_scaled)

    dataG = TGraphErrors(data_x.shape[0], data_x, data_y, data_xerr, data_yerr)
    dataG.SetMarkerStyle(21)
    dataG.GetXaxis().SetRangeUser(-3.0, 4.0)
    dataG.GetYaxis().SetRangeUser(0, 60)
    dataG.GetXaxis().SetTitle('Q-Value [MeV]')
    dataG.GetXaxis().CenterTitle()
    dataG.GetYaxis().SetTitle('Counts')
    dataG.GetYaxis().CenterTitle()
    dataG.SetTitle(r'')
    dataG.Draw('AP')

    simG = TGraph(sim_x.shape[0], sim_x, sim_y_scaled)
    simG.SetMarkerStyle(8)
    simG.SetMarkerColor(2)
    simG.GetXaxis().SetRangeUser(-3.0, 4.0)
    simG.GetYaxis().SetRangeUser(0, 60)
    simG.GetXaxis().SetTitle('Q-Value [MeV]')
    simG.GetXaxis().CenterTitle()
    simG.GetYaxis().SetTitle('Counts')
    simG.GetYaxis().CenterTitle()
    simG.SetTitle(r'')
    simG.Draw('psame')

    text = TLatex(2, 35, 'MSE = ' + '{0:.2f}'.format(mse));
    text.SetTextColor(1)
    text.SetTextSize(0.04)
    text.Draw()

    input_tfile.Close()

    # Output Files
    output_root = TFile(output_file, "recreate")

    c1.Write()

    output_root.Close()

    return mse

if __name__ == '__main__':
    print(analysis(sys.argv[1], sys.argv[2]))
