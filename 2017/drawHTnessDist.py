# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
import sys


def main():

	topness_hh,topness_tt = [],[]
	higgsness_hh,higgsness_tt = [],[]

	hh = open('test_signal.csv','r')
	hhline = csv.reader(hh)

	tt = open('test_background.csv','r')
	ttline = csv.reader(tt)


	error_value = 0
	line_num = 0
	for line in hhline:
		if 'Bjet1_pT' in line: continue
#		if (float(line[-5]) == 0 or float(line[-6]) == 0): continue
		topness_hh.append(np.log10(float(line[-5])))
		higgsness_hh.append(np.log10(float(line[-6])))
		if (np.log10(float(line[-6])) > 0): line_num += 1

#	error_value = 0
#	line_num = 0
	for line in ttline:
		if 'Bjet1_pT' in line: continue
#		if (float(line[-5]) == 0 or float(line[-6]) == 0): continue
#		line_num += 1
#		print(np.log10(float(line[-2])))
#		topness_tt.append(np.log10((line[-2])))
#		higgsness_tt.append(np.log10((line[-3])))
#		try:

#		print("higgsness at tt bar : ", line[-3])
		topness_tt.append(np.log10(float(line[-5])))
		higgsness_tt.append(np.log10(float(line[-6])))
#		if (np.log10(float(line[-3])) > 0): line_num += 1
#		except ValueError:
#			error_value += 1


	c = TCanvas("c","Histogram in (H,T) Space TT Events_SL", 900,600)

	h2 = TH2D("h2","Histogram in (H,T) Space HH and tt Events at Single leptonic decay modes", 100,0,12, 50,0,6)

	for i in range(len(topness_hh)):
		h2.Fill(higgsness_hh[i],topness_hh[i])
	for i in range(len(topness_hh)):
		h2.Fill(higgsness_tt[i],topness_tt[i])

	c.cd()

	h2.GetXaxis().SetTitle("Log H")
	h2.GetYaxis().SetTitle("Log T")
	h2.Draw("colz")

	c.SaveAs("htness_tt.png")

	print("############################")
	print("Value Errors: %i"%(error_value))
	print("Number of Events: %i"%(line_num))
	print("############################")
	

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()


