# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
import sys


def main():

#	norm_hh,prev_hh = [],[]
        norm_hh, norm_tt = [],[]
	norm_jj = []
	norm_kk = []
	norm_rr = []

	nn = open('test_background.csv','r')
	nnline = csv.reader(nn)

#	pp = open('error_hh.csv','r')
#	ppline = csv.reader(pp)

        tt = open('test_background.csv','r')
        ttline = csv.reader(tt)

	jj = open('test_background.csv', 'r')
	jjline = csv.reader(jj)

	kk = open('test_background.csv', 'r')
	kkline = csv.reader(kk)

	rr = open('test_background.csv', 'r')
	rrline = csv.reader(rr)

	for line in nnline:
		if 'Bjet1_pT' in line: continue
		norm_hh.append(float(line[-13])) # -2:topness, -3:higgsness

#	for line in nnline:
#		if 'jet1_mass' in line: continue
#		norm_hh.append(np.log10(float(line[-2])))

#        for line in ttline:
#                if 'jet1_mass' in line: continue
#                norm_tt.append(np.log10(float(line[-2])))


        for line in ttline:
                if 'Bjet1_pT' in line: continue
                norm_tt.append(float(line[-9]))
               
	for line in jjline:
		if 'Bjet1_pT' in line: continue
		norm_jj.append(float(line[-2]))

	for line in kkline:
		if 'Bjet1_pT' in line: continue
		norm_kk.append(float(line[-30]))
		
	for line in rrline:
		if 'Bjet1_pT' in line: continue
		norm_rr.append(float(line[-19]))


	c = TCanvas("c","1D Distribution on Topness", 900,600)
#	c = TCanvas("c","1D Distribution on 2W Mass", 900,600)

#	h1 = TH1D("h1","1D Distribution on Topness_hh", 50,0,5)
#	h2 = TH1D("h2","1D Distribution on Topness_hh", 50,0,5)

#	h1 = TH1D("h1","1D Distribution on Higgsness", 50,0,7)
#	h2 = TH1D("h2","1D Distribution on Higgsness", 50,0,7)

#	h1 = TH1D("h1", "1D Distribution on two W Mass", 50, 0, 500)
#	h2 = TH1D("h2", "1D Distribution on two W  Mass", 50, 0, 500)

#	h1 = TH1D("h1", "1D Distribution on W Mass (lep+nu)", 30, 0, 300)
#	h2 = TH1D("h2", "1D Distribution on W Mass (lep+nu)", 30, 0, 300)

#	h1 = TH1D("h1", "1D Distribution on W Mass (two jets)", 30, 0, 300)
#	h2 = TH1D("h2", "1D Distribution on W Mass (two jets)", 30, 0, 300)

#	h1 = TH1D("h1", "1D Distribution on Expectation Eta ", 140, -7, 7)
#	h2 = TH1D("h2", "1D Distribution on Expectation Eta ", 140, -7, 7)

#	h1 = TH1D("h1", "1D Distributions on Higgs Mass", 30, 0, 300)
#	h2 = TH1D("h2", "1D Distributions on Higgs Mass", 30, 0, 300)
#	h3 = TH1D("h3", "1D Distributions on Higgs Mass", 30, 0, 300)
#
#	h4 = TH1D("h4", "1D Distributions on Higgs Mass", 30, 0, 300)
#	h5 = TH1D("h5", "1D Distributions on Higgs Mass", 30, 0, 300)
#	h6 = TH1D("h6", "1D Distributions on Higgs Mass", 30, 0, 300)

	h1 = TH1D("h1", "1D Distributions on W Mass", 40, 0, 400)
	h2 = TH1D("h2", "1D Distributions on W Mass", 40, 0, 400)
	h3 = TH1D("h3", "1D Distributions on W Mass", 40, 0, 400)	

	h4 = TH1D("h4", "1D Distribution on W Mass", 40, 0, 400)
	h5 = TH1D("h5", "1D Distribution on W Mass", 40, 0, 400)
	
	h6 = TH1D("h6", "1D Distribution on W Mass", 40, 0, 400)

#        h1 = TH1D("h1","1D Distribution on Delta R of lep1 and jet1", 50,0,10)
#        h2 = TH1D("h2","1D Distribution on Delta R of lep1 and jet1", 50,0,10)

#	h1 = TH1D("h1", "Number of Jets", 15, 0, 15)
#	h2 = TH1D("h2", "Number of Jets", 15, 0, 15)
#	h3 = TH1D("h3", "Number of Jets", 15, 0, 15)

#	for i in range(len(prev_hh)):
#		h1.Fill(norm_hh[i])
#		h2.Fill(prev_hh[i])

        for i in range(len(norm_hh)):
                h1.Fill(norm_hh[i])
                h2.Fill(norm_tt[i])
		h3.Fill(norm_jj[i])
		h4.Fill(norm_kk[i])
		h5.Fill(norm_rr[i])

#		print("eta value : ", norm_hh[i])

	c.cd()

#	h6.Add(h4, 1)
	h6.Add(h3, 1)
	h6.Add(h2, 1)
	h6.Add(h1, 1)

	

#	h2.GetXaxis().SetTitle("Log(Topness)")
#	h2.GetYaxis().SetTitle("Entry")
#	h1.SetLineColor(kBlack)
#	h2.SetLineColor(kGray)

#	lst_hh = TLegend(0.80,0.65,0.95,0.75)
#	lst_hh.AddEntry(h1,"norm","l")
#	lst_hh.AddEntry(h2,"prev","l")

#	h2.Draw()
#	h1.Draw("same")
#	lst_hh.Draw()

#-----------------------------------------------------------------------------------------------
#	h5.GetXaxis().SetTitle("Mass (GeV)")
#	h5.GetYaxis().SetTitle("Entry")
#	h5.SetLineColor(kViolet)
#	h6.SetLineColor(kBlue)
#	h5.SetLineWidth(3)
#	h6.SetLineWidth(3)
#
#	lst_hh = TLegend(0.80,0.65,0.95,0.75)
#	lst_hh.AddEntry(h5, "Resolved", "l")
#	lst_hh.AddEntry(h6, "Boosted", "l")
#
#	h5.Draw()
#	h6.Draw("same")
#	lst_hh.Draw()
#-----------------------------------------------------------------------------------------------
	h1.GetXaxis().SetTitle("Mass (GeV)")
        h1.GetYaxis().SetTitle("Entry")
        h1.SetLineColor(kBlue)
        h2.SetLineColor(kRed)
        h1.SetLineWidth(3)
        h2.SetLineWidth(3)

	h2.GetXaxis().SetTitle("Mass (GeV)")
	h2.GetYaxis().SetTitle("Entry")

	h3.GetXaxis().SetTitle("Mass (GeV)")
	h3.GetYaxis().SetTitle("Entry")
	h3.SetLineColor(kOrange)
	h3.SetLineWidth(3)

#        h4.GetXaxis().SetTitle("Mass (GeV)")
#        h4.GetYaxis().SetTitle("Entry")
#        h4.SetLineColor(kViolet)
#        h4.SetLineWidth(3)


        lst_hh = TLegend(0.80,0.65,0.95,0.75)

#	lst_hh.AddEntry(h1, "H to bb", "l")
#	lst_hh.AddEntry(h2, "H to lvqq (MET used)", "l")

	lst_hh.AddEntry(h1, "Off-shell W", "l")
	lst_hh.AddEntry(h2, "On-shell W", "l")

#	lst_hih.AddEntry(h3, "H to lvqq (Neutrino opt used)","l")        

#        lst_hh.AddEntry(h1, "Hadronic", "l")
#        lst_hh.AddEntry(h2, "Leptonic", "l")


#	lst_hh.AddEntry(h1, "t to blv (MET used)", "l")
#	lst_hh.AddEntry(h2, "t to bqq", "l")
#	lst_hh.AddEntry(h3, "t to blv (Neutrino opt used)","l")
#	lst_hh.AddEntry(h4, "FatJet")

        h3.Draw()
	h1.Draw("same")
	h2.Draw("same")
#	h1.Draw("same")
        lst_hh.Draw()



	c.SaveAs("230815_Higgsness.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()



