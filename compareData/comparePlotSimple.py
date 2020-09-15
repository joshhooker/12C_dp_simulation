import numpy as np

from ROOT import gROOT, TCanvas, TColor, TFile, TGaxis, TGraph, TGraphErrors, TPad
from ROOT import kBlack, kBlue, kRed
from array import array

c1 = TCanvas('c1', 'Compare Simulation to Real Data', 200, 10, 700, 500)

data_x, data_y, data_yerr = np.loadtxt("CQvalueTable.txt", skiprows=1, unpack=True)
data_x = data_x.flatten('C'); data_x = data_x.astype('float')
data_y = data_y.flatten('C'); data_y = data_y.astype('float')
data_xerr = np.zeros(data_x.shape[0])
data_yerr = data_yerr.flatten('C'); data_yerr = data_yerr.astype('float')

jeromeGR = TGraphErrors(data_x.shape[0], data_x, data_y, data_xerr, data_yerr)
jeromeGR.SetMarkerStyle(8)
jeromeGR.GetXaxis().SetRangeUser(-6.5, 4.0)
# jeromeGR.GetYaxis().SetRangeUser(0, 40)
jeromeGR.GetXaxis().SetTitle('Q-Value [MeV]')
jeromeGR.GetXaxis().CenterTitle()
jeromeGR.GetYaxis().SetTitle('Counts')
jeromeGR.GetYaxis().CenterTitle()
jeromeGR.SetTitle('')
jeromeGR.Draw('AP')

inFile = TFile.Open("sim.root", "READ")
calcQValueHist = inFile.Get("calcQValue")
i_size = calcQValueHist.GetSize()
xaxis = calcQValueHist.GetXaxis()
scale = 1
highBinContent = 0
sim_x = []
sim_y = []
for i in range(i_size):
    binCenter = xaxis.GetBinCenter(i)
    binContent = calcQValueHist.GetBinContent(i)
    sim_x.append(binCenter)
    sim_y.append(binContent)
    if binContent > highBinContent:
        highBinContent = binContent

scale = 50/highBinContent

sim_x = np.array(sim_x)
sim_y = np.array(sim_y)
sim_y *= scale
simGR = TGraph(sim_x.shape[0], sim_x, sim_y)
simGR.SetLineColor(kRed)
simGR.Draw("LSAME")


# calcQValueHist.Scale(scale)
# calcQValueHist.SetLineColor(kRed)
# calcQValueHist.Draw('LSAME')

# TCanvas.Update() draws the frame, after which one can change it
c1.Update()
c1.GetFrame().SetBorderSize(12)
c1.Modified()
c1.Update()

# Output Files
output_root = TFile("output.root", "recreate")

c1.Write()

output_root.Close()