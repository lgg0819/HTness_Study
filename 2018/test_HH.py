# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
import sys

import pandas as pd
import datetime
import ast
import re

from scipy.optimize import minimize, NonlinearConstraint


def main():

	h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876	

	sig_h = 51.18
	sig_bb = 33.8
	sig_on,sig_off = 40.19,26.86
	peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )
	


	list_FpT, list_FEta, list_FPhi, list_FMass = [], [], [], [] # for FatJet candidates, p1
	list_bpT, list_bEta, list_bPhi, list_bMass = [], [], [], [] # for bJet candidates, p3
	list_jpT, list_jEta, list_jPhi, list_jMass = [], [], [], [] # for Jet candidates, p4
	list_lpT, list_lEta, list_lPhi, list_lMass = [], [], [], [] # for lepton, p5
	list_MpT, list_MEta, list_MPhi, list_MMass = [], [], [], [] # for MET, p6
	list_event = []
	list_case = []

	data = []
	file_name = 'HH_Object_Selection.csv'

	# 1_ = FatJet    2_ = Event case    3_ = BJet candidates    4_ = Jet candidates    5_ = lepton    6_ = MET    7 = Event number
	# _1 = pT,   _2 = Eta,   _3 = Phi,   _4 = Mass


	p1_1 = open(file_name, 'r')
	p1_2 = open(file_name, 'r')
	p1_3 = open(file_name, 'r')
	p1_4 = open(file_name, 'r')

	p2 = open(file_name, 'r')

	p3_1 = open(file_name, 'r')
	p3_2 = open(file_name, 'r')
	p3_3 = open(file_name, 'r')
	p3_4 = open(file_name, 'r')

	p4_1 = open(file_name, 'r')
	p4_2 = open(file_name, 'r')
	p4_3 = open(file_name, 'r')
	p4_4 = open(file_name, 'r')

	p5_1 = open(file_name, 'r')
	p5_2 = open(file_name, 'r')
	p5_3 = open(file_name, 'r')
	p5_4 = open(file_name, 'r')

	p6_1 = open(file_name, 'r')
	p6_2 = open(file_name, 'r')
	p6_3 = open(file_name, 'r')
	p6_4 = open(file_name, 'r')

	p7 = open(file_name, 'r')

	p1_1line = csv.reader(p1_1)
	p1_2line = csv.reader(p1_2)
	p1_3line = csv.reader(p1_3)
	p1_4line = csv.reader(p1_4)

	p2line = csv.reader(p2)

	p3_1line = csv.reader(p3_1)
	p3_2line = csv.reader(p3_2)
	p3_3line = csv.reader(p3_3)
	p3_4line = csv.reader(p3_4)

	p4_1line = csv.reader(p4_1)
	p4_2line = csv.reader(p4_2)
	p4_3line = csv.reader(p4_3)
	p4_4line = csv.reader(p4_4)

	p5_1line = csv.reader(p5_1)
	p5_2line = csv.reader(p5_2)
	p5_3line = csv.reader(p5_3)
	p5_4line = csv.reader(p5_4)

	p6_1line = csv.reader(p6_1)
	p6_2line = csv.reader(p6_2)
	p6_3line = csv.reader(p6_3)
	p6_4line = csv.reader(p6_4)

	p7line = csv.reader(p7)


#	df_divide = pd.read_csv(file_name)
#	divide_columns = ['jet_cand_pt', 'jet_cand_eta', 'jet_cand_phi', 'jet_cand_mass']
#	for col in divide_columns:
#		df_divide[col] = df_divide[col].apply(ast.literal_eval)

#	for row in df_divide.iteruples(index = True):
#		pt_list = row.jet_cand_pt
#		eta_list = row.jet_cand_eta
#		phi_list = row.jet_cand_phi
#		mass_list = row.jet_cand_mass



	for line in p1_1line:
		if 'fatjet_pt' in line: continue
		list_FpT.append(line[-22])

	for line in p1_2line:
		if 'fatjet_pt' in line: continue
		list_FEta.append(line[-21])

	for line in p1_3line:
		if 'fatjet_pt' in line: continue
		list_FPhi.append(line[-20])

	for line in p1_4line:
		if 'fatjet_pt' in line: continue
		list_FMass.append(line[-19])



	for line in p2line:
		if 'fatjet_pt' in line: continue
		list_case.append(float(line[-18]))




	for line in p3_1line:
		if 'fatjet_pt' in line: continue
		list_bpT.append(line[-17])

	for line in p3_2line:
		if 'fatjet_pt' in line: continue
		list_bEta.append(line[-16])

	for line in p3_3line:
		if 'fatjet_pt' in line: continue
		list_bPhi.append(line[-15])

	for line in p3_4line:
		if 'fatjet_pt' in line: continue
		list_bMass.append(line[-14])


	for line in p4_1line:
		if 'fatjet_pt' in line: continue
#		list_jpT.append(float(line[-13]))
		list_jpT.append(line[-13])

	for line in p4_2line:
		if 'fatjet_pt' in line: continue
#		list_jEta.append(float(line[-12]))
		list_jEta.append(line[-12])

	for line in p4_3line:
		if 'fatjet_pt' in line: continue
#		list_jPhi.append(float(line[-11]))
		list_jPhi.append(line[-11])

	for line in p4_4line:
		if 'fatjet_pt' in line: continue
#		list_jMass.append(float(line[-10]))
		list_jMass.append(line[-10])


	for line in p5_1line:
		if 'fatjet_pt' in line: continue
		list_lpT.append(float(line[-9]))

	for line in p5_2line:
		if 'fatjet_pt' in line: continue
		list_lEta.append(float(line[-8]))

	for line in p5_3line:
		if 'fatjet_pt' in line: continue
		list_lPhi.append(float(line[-7]))

	for line in p5_4line:
		if 'fatjet_pt' in line: continue
		list_lMass.append(float(line[-6]))


	for line in p6_1line:
		if 'fatjet_pt' in line: continue
		list_MpT.append(float(line[-5]))

	for line in p6_2line:
		if 'fatjet_pt' in line: continue
		list_MEta.append(float(line[-4]))

	for line in p6_3line:
		if 'fatjet_pt' in line: continue
		list_MPhi.append(float(line[-3]))

	for line in p6_4line:
		if 'fatjet_pt' in line: continue
		list_MMass.append(float(line[-2]))


	for line in p7line:
		if 'fatjet_pt' in line: continue
		list_event.append(float(line[-1]))


	passed_event = 0
	for i in range(len(list_event)):

		if (list_case[0] == 1): continue
#       df_divide = pd.read_csv(file_name)
#       divide_columns = ['jet_cand_pt', 'jet_cand_eta', 'jet_cand_phi', 'jet_cand_mass']
#       for col in divide_columns:
#               df_divide[col] = df_divide[col].apply(ast.literal_eval)

#       for row in df_divide.iteruples(index = True):
#               pt_list = row.jet_cand_pt
#               eta_list = row.jet_cand_eta
#               phi_list = row.jet_cand_phi
#               mass_list = row.jet_cand_mass

#		for j in list_jpT:
		df_pT = ast.literal_eval(list_jpT[i])
		df_Eta = ast.literal_eval(list_jEta[i])
		df_Phi = ast.literal_eval(list_jPhi[i])
		df_Mass = ast.literal_eval(list_jMass[i])

		df_bpT = ast.literal_eval(list_bpT[i])
		df_bEta = ast.literal_eval(list_bEta[i])
		df_bPhi = ast.literal_eval(list_bPhi[i])
		df_bMass = ast.literal_eval(list_bMass[i])


#		print(df_pT)

		Jet1, Jet2, BJet1, BJet2 = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()
		lep_p4, met_p4 = TLorentzVector(), TLorentzVector()
		H1, H2, t1, t2 = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()

#		BJet1.SetPtEtaPhiM(list_b1pT[i], list_b1Eta[i], list_b1Phi[i], list_b1Mass[i])
#		BJet2.SetPtEtaPhiM(list_b2pT[i], list_b2Eta[i], list_b2Phi[i], list_b2Mass[i])

		lep_p4.SetPtEtaPhiM(list_lpT[i], list_lEta[i], list_lPhi[i], list_lMass[i])
		met_p4.SetPtEtaPhiM(list_MpT[i], list_MEta[i], list_MPhi[i], list_MMass[i])
		W_lep = lep_p4 + met_p4

		TCand1, TCand2 = TLorentzVector(), TLorentzVector()
		cand_pt, cand_eta, cand_phi, cand_mass = [], [], [], []
		cand_dM1_1, cand_dM1_2, cand_dM1_3 = [], [], []
		cand_dM2_1, cand_dM2_2, cand_dM2_3 = [], [], []
		total_dM = []
		dR_jj = []
		had_pt1, had_eta1, had_phi1, had_mass1 = [], [], [], []
		had_pt2, had_eta2, had_phi2, had_mass2 = [], [], [], []
		Propagation_error = []

		if (len(df_bpT) > 1):
			BJet1.SetPtEtaPhiM(df_bpT[0], df_bEta[0], df_bPhi[0], df_bMass[0])
			BJet2.SetPtEtaPhiM(df_bpT[1], df_bEta[1], df_bPhi[1], df_bMass[1])
		else: continue

#		if (len(df_bpT) == 1):
#			BJet1.SetPtEtaPhiM(df_bpT[0], df_bEta[0], df_bPhi[0], df_bMass[0])
#			BJet2.SetPtEtaPhiM(0, 0, 0, 0)
#		if (len(df_bpT) > 1):
#			BJet1.SetPtEtaPhiM(df_bpT[1], df_bEta[1], df_bPhi[1], df_bMass[1])
#			for k in range(len(df_bpT)):
#				if (max(df_bMass) == df_bMass[k]):
#					BJet1.SetPtEtaPhiM(df_bpT[k], df_bEta[k], df_bPhi[k], df_bMass[k])
#					BJet2.SetPtEtaPhiM(0, 0, 0, 0)

		if (len(df_pT) == 1):
#			if (df_Mass[0] < w_mass - 15): continue
			Jet1.SetPtEtaPhiM(df_pT[0], df_Eta[0], df_Phi[0], df_Mass[0])
			Jet2.SetPtEtaPhiM(0, 0, 0, 0)
			if (abs(df_Mass[0] - z_mass) < 10): continue
#		if (len(df_pT) >= 2):
#			if (df_Mass[0] < peak_off + 15 and df_Mass[1] < peak_off + 15):
#			Jet1.SetPtEtaPhiM(df_pT[0], df_Eta[0], df_Phi[0], df_Mass[0])
#			Jet2.SetPtEtaPhiM(df_pT[1], df_Eta[1], df_Phi[1], df_Mass[1])
#			if (abs((Jet1 + Jet2).M() - z_mass) < 10): continue
#
#
#			if (df_Mass[0] >= peak_off + 15 or df_Mass[1] >= peak_off + 15):
#				total_dM = [abs(w_mass - df_Mass[0]) + abs(peak_off - W_lep.M()), abs(peak_off - df_Mass[0]) + abs(w_mass - W_lep.M()), abs(w_mass - df_Mass[1]) + abs(peak_off - W_lep.M()), abs(peak_off - df_Mass[1]) + abs(w_mass - W_lep.M()) ]
#				if (min(total_dM) == total_dM[0] or min(total_dM) == total_dM[1]):
#					Jet1.SetPtEtaPhiM(df_pT[0], df_Eta[0], df_Phi[0], df_Mass[0])
#					Jet2.SetPtEtaPhiM(0, 0, 0, 0)
#				if (min(total_dM) == total_dM[2] or min(total_dM) == total_dM[3]):
#					Jet1.SetPtEtaPhiM(0, 0, 0, 0)
#					Jet2.SetPtEtaPhiM(df_pT[1], df_Eta[1], df_Phi[1], df_Mass[1])



#			elif (df_Mass[0] >= peak_off + 15 and df_Mass[1] < peak_off + 15):
#				Jet1.SetPtEtaPhiM(df_pT[0], df_Eta[0], df_Phi[0], df_Mass[0])
#				Jet2.SetPtEtaPhiM(0, 0, 0, 0)
#				if (abs(df_Mass[0] - z_mass) < 10): continue
#				if (W_lep.M() < peak_off + 15): continue
#
#			elif (df_Mass[0] < peak_off + 15 and df_Mass[1] >= peak_off + 15):
#				Jet1.SetPtEtaPhiM(0, 0, 0, 0)
#				Jet2.SetPtEtaPhiM(df_pT[1], df_Eta[1], df_Phi[1], df_Mass[1])
#				if (abs(df_Mass[1] - z_mass) < 10): continue
#				if (W_lep.M() < peak_off + 15): continue
#			else:
#				continue
		if (len(df_pT) >= 2):
			led_onejet = 0
			for j in range(len(df_pT)):
				E_Jet1 = 0
	
				if (j >= len(df_pT) - 1): continue
				if (df_Mass[j] > w_mass - 15):
#					led_onejet = 1
#					Jet1.SetPtEtaPhiM(df_pT[j], df_Eta[j], df_Phi[j], df_Mass[j])
#					Jet2.SetPtEtaPhiM(0, 0, 0, 0)
					TCand1.SetPtEtaPhiM(df_pT[j], df_Eta[j], df_Phi[j], df_Mass[j])
					TCand2.SetPtEtaPhiM(0, 0, 0, 0)
					if (abs(z_mass - (TCand1 + TCand2).M() ) <= 10): continue
#					if (abs( (TCand1 + TCand2).M() - W_lep.M()) <= 15): continue
					cand_pt.append( (TCand1 + TCand2).Pt() )
					cand_eta.append( (TCand1 + TCand2).Eta() )
					cand_phi.append( (TCand1 + TCand2).Phi() )
					cand_mass.append( (TCand1 + TCand2).M() )

					W_cand = TCand1 + TCand2
#					dR_jj.append(W_cand.DeltaR(W_lep))
					dR_jj.append(TCand1.DeltaR(TCand2))

					had_pt1.append(TCand1.Pt())
					had_eta1.append(TCand1.Eta())
					had_phi1.append(TCand1.Phi())
					had_mass1.append(TCand1.M())

					had_pt2.append(0)
					had_eta2.append(0)
					had_phi2.append(0)
					had_mass2.append(0)

					continue

#				if (led_onejet == 1): break

				TCand1.SetPtEtaPhiM(df_pT[j], df_Eta[j], df_Phi[j], df_Mass[j])
				E_Jet1 = TCand1.Pt() * math.cosh(TCand1.Eta())
				for k in range(len(df_pT)):
					E_Jet2 = 0
					if (k <= j): continue
					if (df_Mass[k] > w_mass - 15): continue
						
					TCand2.SetPtEtaPhiM(df_pT[k], df_Eta[k], df_Phi[k], df_Mass[k])
					E_Jet2 = TCand2.Pt() * math.cosh(TCand2.Eta())
					if (abs(z_mass - (TCand1 + TCand2).M() ) <= 10): continue
#					if (abs( (TCand1 + TCand2).M() - W_lep.M()) <= 15): continue
					cand_pt.append( (TCand1 + TCand2).Pt() )
					cand_eta.append( (TCand1 + TCand2).Eta() )
					cand_phi.append( (TCand1 + TCand2).Phi() )
					cand_mass.append( (TCand1 + TCand2).M() )

					cand_dM1_1.append( abs( w_mass - (TCand1 + TCand2).M() ) + abs( peak_off - W_lep.M() ) ) 
					cand_dM1_2.append( abs( peak_off - (TCand1 + TCand2).M() ) + abs( w_mass - W_lep.M() ) )
					cand_dM1_3.append( abs( w_mass - (TCand1 + TCand2).M() ) + abs( w_mass - W_lep.M() ) )

					W_cand = TCand1 + TCand2
#					if (TCand1.DeltaR(TCand2) <= 0.4): continue
	
					dR_jj.append(TCand1.DeltaR(TCand2))
#					dR_jj.append(W_cand.DeltaR(W_lep))

					had_pt1.append(TCand1.Pt())
					had_eta1.append(TCand1.Eta())
					had_phi1.append(TCand1.Phi())
					had_mass1.append(TCand1.M())
	
					had_pt2.append(TCand2.Pt())
					had_eta2.append(TCand2.Eta())
					had_phi2.append(TCand2.Phi())
					had_mass2.append(TCand2.M())



			for j in range(len(dR_jj)):
				if (min(dR_jj) == dR_jj[j]):
					Jet1.SetPtEtaPhiM(had_pt1[j], had_eta1[j], had_phi1[j], had_mass1[j])
					Jet2.SetPtEtaPhiM(had_pt2[j], had_eta2[j], had_phi2[j], had_mass2[j])

#			for j in range(len(df_pT)):
				

#			if (len(cand_dM1_1) == 0 or len(cand_dM1_2) == 0 or len(cand_dM1_3) == 0): continue
#			total_dM = [min(cand_dM1_1), min(cand_dM1_2), min(cand_dM1_3)]
#
#			if (min(total_dM) == min(cand_dM1_1)):
#				for j in range(len(cand_dM1_1)):
#					if (min(cand_dM1_1) == cand_dM1_1[j]):
#						Jet1.SetPtEtaPhiM(had_pt1[j], had_eta1[j], had_phi1[j], had_mass1[j])
#						Jet2.SetPtEtaPhiM(had_pt2[j], had_eta2[j], had_phi2[j], had_mass2[j])
#						break
#
#			if (min(total_dM) == min(cand_dM1_2)):
#				for j in range(len(cand_dM1_2)):
#					if (min(cand_dM1_2) == cand_dM1_2[j]):
#						Jet1.SetPtEtaPhiM(had_pt1[j], had_eta1[j], had_phi1[j], had_mass1[j])
#						Jet2.SetPtEtaPhiM(had_pt2[j], had_eta2[j], had_phi2[j], had_mass2[j])
#						break
#
#			if (min(total_dM) == min(cand_dM1_3)):
#				for j in range(len(cand_dM1_3)):
#					if (min(cand_dM1_3) == cand_dM1_3[j]):
#						Jet1.SetPtEtaPhiM(had_pt1[j], had_eta1[j], had_phi1[j], had_mass1[j])
#						Jet2.SetPtEtaPhiM(had_pt2[j], had_eta2[j], had_phi2[j], had_mass2[j])
#						break

		W_had = Jet1 + Jet2


#		if (len(list_jpT) == -1):
#		Jet1.SetPtEtaPhiM(df_pT[0], df_Eta[0], df_Phi[0], df_Mass[0])
#		Jet2.SetPtEtaPhiM(0, 0, 0, 0)

#		print("Jet1 pT : {}    Jet1 Eta : {}    Jet1 Phi : {}    Jet1 Mass : {} ".format(Jet1.Pt(), Jet1.Eta(), Jet1.Phi(), Jet1.M()))
#		if (len(list_jpT) < -1):
#			Jet1.SetPtEtaPhiM(float(list_jpT[-1]), float(list_jEta[-1]), float(list_jPhi[-1]), float(list_jMass[-1]))
#			Jet2.SetPtEtaPhiM(float(list_jpT[-2]), float(list_jEta[-2]), float(list_jPhi[-2]), float(list_jMass[-2]))





		W_on, W_off = TLorentzVector(), TLorentzVector()

#		if (Jet2.Pt() == 0): continue

#		if (BJet1.DeltaR(lep_p4) < 0.4 and BJet2.DeltaR(lep_p4) < 0.4): continue

#		if (Jet1.DeltaR(BJet1) < 2.0 or Jet1.DeltaR(BJet2) < 2.0): continue
#		if (Jet2.DeltaR(BJet1) < 2.0 or Jet2.DeltaR(BJet2) < 2.0): continue
#		if (Jet1.DeltaR(Jet2) > 1.5): continue		


#		if (W_lep.DeltaR(W_had) < 0.8 and Jet1.DeltaR(Jet2) > 0.4):
#		if (W_lep.DeltaR(W_had) < 0.8):
#			H1 = BJet1 + BJet2
#			H2 = W_lep + W_had
	
#			if (W_lep.DeltaR(W_had) < 0.8 and Jet2.Pt() == 0):
#					H1 = BJet1 + BJet2
#					H2 = W_lep + W_had



		if (Jet1.Pt() == 0 or Jet2.Pt() == 0):
			if (abs(W_had.M() - w_mass) < 10):
				t1 = BJet1 + W_lep
				t2 = BJet2 + W_had
		else:
			t1 = BJet2 + W_lep
			t2 = BJet1 + W_had

		Chi_W = []
		Chi_W.append((W_had.M()-w_mass)**2/sig_on**2 + (W_lep.M()-peak_off)**2/sig_off**2)
		Chi_W.append((W_lep.M()-w_mass)**2/sig_on**2 + (W_had.M()-peak_off)**2/sig_off**2)


#		if (min(Chi_W) == Chi_W[0]):
		if (peak_off + 15 >= W_lep.M()):
#		if (abs(peak_off - W_lep.M()) <= 10):
			W_on = W_had
			W_off = W_lep

#		if (min(Chi_W) == Chi_W[1]):
		if (peak_off + 15 < W_lep.M()):
#		if (abs(peak_off - W_had.M()) <= 10):
			W_on = W_lep
			W_off = W_had


#		if (W_off.Pt() > 200): continue
#		if (W_lep.DeltaR(W_had) < 0.8):

		if (W_had.M() < w_mass + 20 and W_lep.M() < w_mass + 20):
			if (abs(W_had.M() - W_lep.M()) > 15):
				H1 = BJet1 + BJet2
				H2 = W_lep + W_had
#			H2 = W_on + W_off


		# test cut

		if (W_had.Pt() == 0 or W_lep.Pt() == 0): continue
		if (H1.Pt() == 0 or H2.Pt() == 0): continue
#		if (W_had.M() < peak_off - 15 and W_lep.M() < peak_off - 15): continue # off-shell + off-shell


		# Common veto cut
#		if (W_had.M() <= peak_off + 10 and W_lep.M() <= peak_off + 10): continue # off-shell + off-shell
#		if (W_had.M() <= w_mass - 15 and W_lep.M() <= w_mass - 15): continue # off-shell + off-shell

		# veto cut for tt
#		if (t1.Pt() == 0 or t2.Pt() == 0): continue
#		if (W_had.M() < w_mass - 15 or W_lep.M() < w_mass - 15): continue # off-shell exists


		nu_p4 = TLorentzVector()
		W_nu = lep_p4 + nu_p4
		H_opt = TLorentzVector()

		def defineHiggsness(x):
			chi_Higgs = []
			nu_p4.SetPtEtaPhiM(met_p4.Pt(), x, met_p4.Phi(), 0)
			W_nu = lep_p4 + nu_p4
			H_opt = W_had + W_nu

			sig_h = 48.96
			sig_bb = 34.78
			sig_on,sig_off = 20.5,23.5
			peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

			chi_Higgs.append( (H_opt.M()**2-h_mass**2)**2/sig_h**4 + (W_had.M()**2-w_mass**2)**2/sig_on**4 + (W_nu.M()**2-peak_off**2)**2/sig_off**4 )
			chi_Higgs.append( (H_opt.M()**2-h_mass**2)**2/sig_h**4 + (W_nu.M()**2-w_mass**2)**2/sig_on**4 + (W_had.M()**2-peak_off**2)**2/sig_off**4 )

			higgsness = min(chi_Higgs)
			return higgsness

		ini = np.array([0])
		bnd = [(-2.4,2.4)]
		con_eq = lambda x: np.sqrt( 2*met_p4.Pt()*lep_p4.Pt()*( np.cosh(x-lep_p4.Eta()) - np.cos(met_p4.Phi()-lep_p4.Phi()) ) )
		cons = NonlinearConstraint(con_eq, lep_p4.M(), 2*w_mass)
		opt1 = minimize(defineHiggsness, ini, constraints=cons, bounds=bnd, method='SLSQP', options={'ftol':5}) #'disp':True})

		Higgsness = opt1.fun
		if (BJet1.DeltaR(lep_p4) > 0.4 and BJet2.DeltaR(lep_p4) > 0.4):
			W_nu = lep_p4 + nu_p4
			H_opt = W_had + W_nu

		chi_t = []
		tW_nu = TLorentzVector()
		t_opt = TLorentzVector()
		B_opt = TLorentzVector()

		def defineTopness():
			tW_nu = lep_p4 + nu_p4

			sig_tlep = 45.1
			sig_thad = 29.7
			sig_had = 14.66
			sig_lep = 39.77

			chi_t.append( ((tW_nu + BJet1).M()**2-t_mass**2)**2/sig_tlep**4 + ((W_had + BJet2).M()**2-t_mass**2)**2/sig_thad**4 + (tW_nu.M()**2-w_mass**2)**2/sig_lep**4 + (W_had.M()**2-w_mass**2)**2/sig_had**4)
			chi_t.append( ((tW_nu + BJet2).M()**2-t_mass**2)**2/sig_tlep**4 + ((W_had + BJet1).M()**2-t_mass**2)**2/sig_thad**4 + (tW_nu.M()**2-w_mass**2)**2/sig_lep**4 + (W_had.M()**2-w_mass**2)**2/sig_had**4)

			topness = min(chi_t)
			return topness

		Topness = defineTopness()

		if (min(chi_t) == chi_t[0]): B_opt = BJet1
		if (min(chi_t) == chi_t[1]): B_opt = BJet2

		if (Jet1.Pt() == 0 or Jet2.Pt() == 0):
			if (abs(W_had.M() - w_mass) < 10):
				tW_nu = lep_p4 + nu_p4
				t_opt = tW_nu + B_opt
		else:
			tW_nu = lep_p4 + nu_p4
			t_opt = tW_nu + B_opt

		if (opt1.success): passed_event += 1
		else: continue


		data.append([BJet1.Pt(), BJet1.Eta(), BJet1.Phi(), BJet1.M(),
				BJet2.Pt(), BJet2.Eta(), BJet2.Phi(), BJet2.M(),
				Jet1.Pt(), Jet1.Eta(), Jet1.Phi(), Jet1.M(),
				Jet2.Pt(), Jet2.Eta(), Jet2.Phi(), Jet2.M(),
				lep_p4.Pt(), lep_p4.Eta(), lep_p4.Phi(), lep_p4.M(),
				met_p4.Pt(), met_p4.Eta(), met_p4.Phi(), met_p4.M(),
				W_on.Pt(), W_on.Eta(), W_on.Phi(), W_on.M(), # lep/on 
				W_off.Pt(), W_off.Eta(), W_off.Phi(), W_off.M(), # had/off
				H1.Pt(), H1.Eta(), H1.Phi(), H1.M(),
				H2.Pt(), H2.Eta(), H2.Phi(), H2.M(),
				t1.Pt(), t1.Eta(), t1.Phi(), t1.M(),
				t2.Pt(), t2.Eta(), t2.Phi(), t2.M(),
				W_lep.DeltaR(W_had), BJet1.DeltaR(BJet2), Higgsness, Topness,
				passed_event, H_opt.M(), t_opt.M(), 0])

	df = pd.DataFrame(data, columns=
				['Bjet1_pT', 'Bjet1_Eta', 'Bjet1_Phi', 'Bjet1_Mass',
				'Bjet2_pT', 'Bjet2_Eta', 'Bjet2_Phi', 'Bjet2_Mass',
				'jet1_pT', 'jet1_Eta', 'jet1_Phi', 'jet1_Mass',
				'jet2_pT', 'jet2_Eta', 'jet2_Phi', 'jet2_Mass',
				'lep_pT', 'lep_Eta', 'lep_Phi', 'lep_Mass',
				'MET_pT', 'MET_Eta', 'MET_Phi', 'MET_Mass',
				'W_lep_pT', 'W_lep_Eta', 'W_lep_Phi', 'W_lep_Mass',
				'W_had_pT', 'W_had_Eta', 'W_had_Phi', 'W_had_Mass',
				'H_bb_pT', 'H_bb_Eta', 'H_bb_Phi', 'H_bb_Mass',
				'H_WW_pT', 'H_WW_Eta', 'H_WW_Phi', 'H_WW_Mass',
				't_lep_pT', 't_lep_Eta', 't_lep_Phi', 't_lep_Mass',
				't_had_pT', 't_had_Eta', 't_had_Phi', 't_had_Mass',
				'dR_WW', 'dR_bb', 'Higgsness', 'Topness',
				'Event_Number', 'H_opt_Mass', 't_opt_Mass', 'target'])

	df.to_csv("test_signal.csv", header=True, index=False)

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()
