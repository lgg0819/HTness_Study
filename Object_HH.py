from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
import datetime
from scipy.optimize import minimize, NonlinearConstraint


def main():

	now1 = datetime.datetime.now()
	print("***********************************************************")
	print("Analysis started time: ", now1.strftime("%Y-%m-%d %H:%M:%S"))
	print("***********************************************************")

	# Signal SL
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2WToLNu2J_node_SM_TuneCUETP8M1_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')
	flist = glob.glob('/u/user/lgg0819/mc/GluGluToHHTo2B2WToLNu2J_node_SM_TuneCUETP8M1_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')

	# TTBar 
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')
	
#	flist = glob.glob('/u/user/lgg0819/mc/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8NANOAODSIMPUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/100000/*')

	previous_hh = 0
	ii,hh,ee = 0,0,0
	passed_event = 0
	h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876

	sig_h = 56.84 # standard deviation for H to WW mass
	sig_bb = 37.47 # standard deviation for H to bb mass
	sig_on,sig_off = 48.09,30.93 # standard deviation for W masses at signal (HH)
	peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) ) # peak value of Off-shell W mass distribution

	sig_t = 44.1 # standard deviation for top mass
	sig_w = 19.1 # standard deviation for W mass at main background (tt)
	sig_thad = 72.68
	sig_tlep = 59.71
	sig_had = 50.25
	sig_lep = 44.86


	hgg,w_on,w_off = [],[],[]
	tt,ww = [],[]
	data = []

	good_passed, overlaps = 0, 0 #for check the b-tag overlappnig exists or not
	overlap1, overlap2, overlap3, overlap_all = 0, 0, 0, 0 #for check each overlap cases
	b_tag1, b_tag2, b_tag3, b_tag4, b_tag_all = 0, 0, 0, 0, 0 #for count the b-tag selections
	cross_over = 0
	count_Cross = 0

	count_btag1, count_btag2, count_btag3, count_btag4 = 0, 0, 0, 0
	count_wsel1, count_wsel2, count_wsel3, count_wsel4 = 0, 0, 0, 0

	count_nFat, count_tFat = 0, 0
	count_dilep, count_lep_veto = 0, 0


	type_error = 0

	zero_range_w = 0
	zero_range_t = 0


	count_case1, count_case2, count_case3, count_case4, count_case5, count_case6 = 0, 0, 0, 0, 0, 0

	count_Won, count_Woff, count_Hww, count_Hbb = 0, 0, 0, 0 # added by 240624
	count_Wlv, count_Wqq, count_tlv, count_tqq = 0, 0, 0, 0 # added by 240624

	check_qq1, check_qq2, check_qq3, check_qq4, check_qq5 = 0, 0, 0, 0, 0 
	check_qq6, check_qq7, check_qq8 = 0, 0, 0 
	check_total = 0


	for fname in flist:


		current_time1 = datetime.datetime.now()
#		if passed_event > 15000: continue
		f = TFile(fname,"read")
		t = f.Get("Events")

		try: event_num = t.GetEntriesFast()
		except AttributeError: 
			print("AttributeError appared")
			continue

		print("#########################################")
		print("Analysis started time: ", now1.strftime("%Y-%m-%d %H:%M:%S"))
		print("Current time: ", current_time1.strftime("%Y-%m-%d %H:%M:%S"), "Number of Passed Events: %i"%(passed_event))
		print("#########################################")

		print("Processing %s..."%(fname))
		ii = ii+1
		print("Have processed %i events..."%(ee))
		print("")
		print("")

		LED_checkpoint = 0

		# Only Mu channel considered at this point
		# HLT conditions
		for i in range(event_num):
			current_time2 = datetime.datetime.now()
			if (passed_event != 0 and passed_event%1500 == 0 and passed_event != LED_checkpoint):
				LED_checkpoint = passed_event
				print("#########################################")
				print("Analysis started time: ", now1.strftime("%Y-%m-%d %H:%M:%S"))
				print("Current time: ", current_time2.strftime("%Y-%m-%d %H:%M:%S"), "Number of Passed Events: %i"%(passed_event))
				print("#########################################")
				print("Have processed %i events..."%(ee))


#			if (hh == 5000): break

			Higgsness, Topness = 0, 0


                        t.GetEntry(i)
			hlt_mu,hlt_el = False,False
			t_hlt_mu, t_hlt_el = False, False

			good_event = False
			lep_veto = False
			led_HH = False ## added by 240530

			# SingleMuon
                        hlt_mu_1 = t.HLT_IsoMu22
                        hlt_mu_2 = t.HLT_IsoTkMu22
                        hlt_mu_3 = t.HLT_IsoMu22_eta2p1
                        hlt_mu_4 = t.HLT_IsoTkMu22_eta2p1
                        hlt_mu_5 = t.HLT_IsoMu24
                        hlt_mu_6 = t.HLT_IsoTkMu24

			# reference of trigger at tt bar : (https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger)

                        if (hlt_mu_1 == 1 or hlt_mu_2 == 1 or hlt_mu_3 == 1 or hlt_mu_4 == 1 or hlt_mu_5 == 1 or hlt_mu_6 == 1): hlt_mu = True
			if (hlt_mu_4 == 1 or hlt_mu_6 == 1): t_hlt_mu = True

			# SingleElectron
                        hlt_el_1 = t.HLT_Ele27_WPTight_Gsf
                        hlt_el_2 = t.HLT_Ele25_eta2p1_WPTight_Gsf
                        hlt_el_3 = t.HLT_Ele27_eta2p1_WPLoose_Gsf

			hlt_el_4 = t.HLT_Ele32_eta2p1_WPTight_Gsf

                        if (hlt_el_1 == 1 or hlt_el_2 == 1 or hlt_el_3 == 1): hlt_el = True
			if (hlt_el_1 == 1 or hlt_el_2 == 1 or hlt_el_4 == 1): t_hlt_el = True

                        flag_good  = t.Flag_goodVertices
                        flag_halo  = t.Flag_globalSuperTightHalo2016Filter
                        flag_hbhen = t.Flag_HBHENoiseFilter
                        flag_hbiso = t.Flag_HBHENoiseIsoFilter
                        flag_dead  = t.Flag_EcalDeadCellTriggerPrimitiveFilter
                        flag_badpf = t.Flag_BadPFMuonFilter
                        flag_ecbad = t.Flag_ecalBadCalibFilter
                        flag_eebad = t.Flag_eeBadScFilter

                        met_filter = False
			tmet_filter = False # MET filter for tt bar at 2016 MC
                        if (flag_good == 1 and flag_halo == 1 and flag_hbhen == 1 and flag_hbiso == 1 and flag_dead == 1 and flag_badpf == 1 and flag_ecbad == 1):
                                met_filter = True
			if (flag_good == 1 and flag_halo == 1 and flag_hbhen == 1 and flag_hbiso == 1 and flag_dead == 1 and flag_badpf == 1):
				tmet_filter = True

			if (hlt_mu == True or hlt_el == True):
				lep_veto = True
				if met_filter: good_event = True

			# =================================================================
			#Remained cut : The lepton isolation, defined as the scalar
			#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
			#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
			# Medium muon discrimination

			# =================================================================

			# Array Definition for FatJet, SubJet, Resolved Jet(thin jet)
			fat_mass, fat_pt, fat_eta, fat_phi = [], [], [], []
			thin_mass, thin_pt, thin_eta, thin_phi = [], [], [], []
			sub_mass, sub_pt, sub_eta, sub_phi = [], [], [], []

			fat_btag, fat_jetid, fat_tau1, fat_tau2 = [], [], [], [] #fat_tau'N': N-Subjettiness
			fat_sdm = [] # SoftDrop mass of FatJet
			thin_btag, thin_jetid, thin_charge = [], [], []
			sub_btag, sub_jetid, sub_charge = [], [], []

			thin_tbtag = [] # added by 240513

			thin_nMu, thin_nEl = [], []



			# =================================================================

			ee += 1
			pt,phi = 0,0

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0
			ljet_mass,ljet_pt,ljet_eta,ljet_phi,ljet_btag = 0,0,0,0,0
			sjet_mass,sjet_pt,sjet_eta,sjet_phi,sjet_btag = 0,0,0,0,0
			tjet_mass,tjet_pt,tjet_eta,tjet_phi = 0,0,0,0
			fjet_mass,fjet_pt,fjet_eta,fjet_phi = 0,0,0,0

			el_mass, el_pt, el_eta, el_phi, el_charge = 0,0,0,0,0

			tmu_mass, tmu_pt, tmu_eta, tmu_phi, tmu_charge = 0,0,0,0,0
			tel_mass, tel_pt, tel_eta, tel_phi, tel_charge = 0,0,0,0,0

			## ----------------------for Muon -------------------------------
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                        temp_charge,temp_medid = [],[]
                        temp_dxy,temp_dz,temp_btag = [],[],[]
                        temp_iso,temp_sip3d,temp_mva = [],[],[]

			## ----------------------for Electron ---------------------------
			elep_mass, elep_pt, elep_eta, elep_phi = [], [], [], []
			elep_charge, elep_medid = [], []
			elep_dxy, elep_dz = [], []
			elep_miniPFRelIso_all, elep_sip3d = [], []
			elep_mvaNoIso_WPL, elep_lostHits = [], []
			## ----------------------for Tau veto ---------------------------
			tau_mass, tau_pt, tau_eta, tau_phi = [], [], [], []
			tau_dz = []			
			tau_jet = [] # Tau_idDeepTau2017v2p1VSjet (added by 240531)
			tau_mu = [] # Tau_idDeepTau2017v2p1VSmu (added by 240829)
			tau_e = [] # Tau_idDeepTau2017v2p1VSe (added by 240829)

			# DeepTau Jet: 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
			# DeepTau Mu: 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight
			# DeepTau e: 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight

			## --------------------- added by 240321 ------------------------
			b_jet_pt1, b_jet_eta1, b_jet_phi1, b_jet_mass1 = 0, 0, 0, 0
			b_jet_pt2, b_jet_eta2, b_jet_phi2, b_jet_mass2 = 0, 0, 0, 0

			ak4_boost_pt1, ak4_boost_eta1, ak4_boost_phi1, ak4_boost_mass1 = 0, 0, 0, 0
			ak4_boost_pt2, ak4_boost_eta2, ak4_boost_phi2, ak4_boost_mass2 = 0, 0, 0, 0

			ak4_res_pt1, ak4_res_eta1, ak4_res_phi1, ak4_res_mass1 = 0, 0, 0, 0
			ak4_res_pt2, ak4_res_eta2, ak4_res_phi2, ak4_res_mass2 = 0, 0, 0, 0


			ak_bres_pt1, ak_bres_eta1, ak_bres_phi1, ak_bres_mass1 = 0, 0, 0, 0
			ak_bres_pt2, ak_bres_eta2, ak_bres_phi2, ak_bres_mass2 = 0, 0, 0, 0
			## --------------------------------------------------------------

                        ##---------------------- added by 240406 ------------------------
                        t_ak4_res_pt1, t_ak4_res_eta1, t_ak4_res_phi1, t_ak4_res_mass1 = 0, 0, 0, 0
                        t_ak4_res_pt2, t_ak4_res_eta2, t_ak4_res_phi2, t_ak4_res_mass2 = 0, 0, 0, 0
                        ## --------------------------------------------------------------


			jet1_pt, jet1_eta, jet1_phi, jet1_mass, jet1_btag = 0, 0, 0, 0, 0
			jet2_pt, jet2_eta, jet2_phi, jet2_mass, jet2_btag = 0, 0, 0, 0, 0
			jet3_pt, jet3_eta, jet3_phi, jet3_mass, jet3_btag = 0, 0, 0, 0, 0
			jet4_pt, jet4_eta, jet4_phi, jet4_mass, jet4_btag = 0, 0, 0, 0, 0

			jet1_id, jet2_id, jet3_id, jet4_id = 0, 0, 0, 0
			jet1_tbtag, jet2_tbtag, jet3_tbtag, jet4_tbtag = 0, 0, 0, 0
			jet1_JEC, jet2_JEC, jet3_JEC, jet4_JEC = 0, 0, 0, 0

			flag_lep,flag_jet1,flag_jet2 = False,False,False   ## for Mu, jet1, jet2
			flag_el = False ## for Electron, add by 240318

			tau_cand_pt, tau_cand_eta, tau_cand_phi, tau_cand_mass = [], [], [], []

                        genweight = t.genWeight
			tightMuon = []
			count_tightM, count_tightE = 0, 0

			green_mu = 0 # green -> green light(meaning passed the HLT)
			for j in range(t.nMuon):
				temp_mass.append(t.Muon_mass.__getitem__(j))
				temp_pt.append(t.Muon_pt.__getitem__(j))
				temp_eta.append(t.Muon_eta.__getitem__(j))
				temp_phi.append(t.Muon_phi.__getitem__(j))
				temp_dxy.append(t.Muon_dxy.__getitem__(j))
				temp_dz.append(t.Muon_dz.__getitem__(j))
				temp_charge.append(t.Muon_charge.__getitem__(j))
				temp_medid.append(t.Muon_mediumId.__getitem__(j))
                                temp_iso.append(t.Muon_miniPFRelIso_all.__getitem__(j))
                                temp_sip3d.append(t.Muon_sip3d.__getitem__(j))
                                temp_mva.append(t.Muon_mvaTTH.__getitem__(j))
				if (t.Muon_tightId.__getitem__(j) == True): 
					tightMuon.append(1)
					count_tightM += 1
				else: tightMuon.append(0)
				if (hlt_mu == True): green_mu = 1  ## for HH
#				if (t_hlt_mu == True): green_mu = 1 ## for tt
			tightEl = []

			green_el = 0 # green -> green light(meaning passed the HLT)
                        for j in range(t.nElectron):
                                elep_mass.append(t.Electron_mass.__getitem__(j))
                                elep_pt.append(t.Electron_pt.__getitem__(j))
                                elep_eta.append(t.Electron_eta.__getitem__(j))
                                elep_phi.append(t.Electron_phi.__getitem__(j))
                                elep_dxy.append(t.Electron_dxy.__getitem__(j))
                                elep_dz.append(t.Electron_dz.__getitem__(j))
                                elep_charge.append(t.Electron_charge.__getitem__(j))
				elep_medid.append(t.Electron_mvaTTH.__getitem__(j))
				elep_miniPFRelIso_all.append(t.Electron_miniPFRelIso_all.__getitem__(j))
				elep_sip3d.append(t.Electron_sip3d.__getitem__(j))
				elep_mvaNoIso_WPL.append(t.Electron_mvaFall17V2noIso.__getitem__(j))
				elep_lostHits.append(t.Electron_lostHits.__getitem__(j))
                                if (t.Electron_mvaTTH > 0.30): 
					tightEl.append(1)
					count_tightE += 1
				else: tightEl.append(0)
				if (hlt_el == True): green_el = 1 ## for HH

#			if (count_tightM > 1 or count_tightE > 1): continue # veto for dileptonic channel

			flag_Mu, flag_El = 0, 0
			flag_tMu, flag_tEl = 0, 0

			l1 = -1
			count_M, count_E = 0, 0
			for k in range(t.nMuon):
#				if (flag_lep == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_medid[k] == True and temp_iso[k] < 0.4 and temp_sip3d[k] < 8 and temp_mva[k] > 0.5 and tightMuon[k] == 1):
					if (hlt_mu == True):

						llep_mass = temp_mass[k]
						llep_pt   = temp_pt[k]
						llep_eta  = temp_eta[k]
						llep_phi  = temp_phi[k]
						llep_charge = temp_charge[k]

						l1 = k
						flag_lep = True
						flag_Mu = 1
						count_M += 1

					if (t_hlt_mu == True):

						tmu_mass = temp_mass[k]
						tmu_pt = temp_pt[k]
						tmu_eta = temp_eta[k]
						tmu_phi = temp_phi[k]
						tmu_charhe = temp_charge[k]

						flag_tMu = 1

			if (count_M > 1): continue

			for k in range(t.nElectron):
#				if (flag_el == True): break
				if (elep_pt[k] > 30 and abs(elep_eta[k]) < 2.4 and abs(elep_dxy[k]) < 0.05 and abs(elep_dz[k]) < 0.1 and elep_medid[k] > 0.30 and elep_miniPFRelIso_all[k] < 0.4 and elep_sip3d[k] <= 8 and elep_mvaNoIso_WPL[k] == True and tightEl[k] == 1):

					if (hlt_el == True):
						el_pt = elep_pt[k]
						el_eta = elep_eta[k]
						el_phi = elep_phi[k]
						el_mass = elep_mass[k]
						el_charge = elep_charge[k]

						flag_el = True
						flag_El = 1
						count_E += 1

					if (t_hlt_el == True):
						tel_pt = elep_pt[k]
						tel_eta = elep_eta[k]
						tel_phi = elep_phi[k]
						tel_mass = elep_mass[k]
						tel_charge = elep_charge[k]

						flag_tEl = 1

			if (count_E > 1): continue

                        flag_tau1, flag_tau2 = False, False
                        for j in range(t.nTau):
                                tau_mass.append(t.Tau_mass.__getitem__(j))
                                tau_pt.append(t.Tau_pt.__getitem__(j))
                                tau_eta.append(t.Tau_eta.__getitem__(j))
                                tau_phi.append(t.Tau_phi.__getitem__(j))
                                tau_dz.append(t.Tau_dz.__getitem__(j))
                                tau_jet.append(t.Tau_idDeepTau2017v2p1VSjet.__getitem__(j))
                                tau_mu.append(t.Tau_idDeepTau2017v2p1VSmu.__getitem__(j))
                                tau_e.append(t.Tau_idDeepTau2017v2p1VSe.__getitem__(j))


                        for j in range(len(tau_pt)):
                                if (tau_pt[j] > 40 and abs(tau_eta[j]) < 2.3 and tau_dz[j] < 0.2):
                                        if (tau_jet[j] >= 16 and tau_mu[j] >= 8 and tau_e[j] >= 32 and hlt_mu == True):
                                                tau_cand_pt.append(tau_pt[j])
                                                tau_cand_eta.append(tau_eta[j])
                                                tau_cand_phi.append(tau_phi[j])
                                                tau_cand_mass.append(tau_mass[j])
			for j in range(len(tau_cand_pt)):
				if (math.sqrt((tau_cand_eta[j] - llep_eta)**2 + (tau_cand_phi[j] - llep_phi)**2) > 0.3): flag_tau1 = True
				if (math.sqrt((tau_cand_eta[j] - el_eta)**2 + (tau_cand_phi[j] - el_phi)**2) > 0.3): flag_tau2 = True

			if (hlt_mu == True and hlt_el == False and flag_tau1 == True): continue
			if (hlt_mu == False and hlt_el == True and flag_tau2 == True): continue 

			# =================================================================
			# Jet discrimination

			met_pt, met_phi = 0, 0
			tmet_pt, tmet_phi = 0, 0

			if (met_filter == True and t.MET_pt > 20):
				met_pt  = t.MET_pt
				met_phi = t.MET_phi



			# =================================================================
			# Preparing for SubJet Selection
			nSub = 0
			flag_sub1 = 0 #for number of SubJets
			flag_sub2 = 0 #for condition that pT > 30 GeV
			for j in range(t.nSubJet):
				sub_mass.append(t.SubJet_mass.__getitem__(j))
				sub_pt.append(t.SubJet_pt.__getitem__(j))
				sub_eta.append(t.SubJet_eta.__getitem__(j))
				sub_phi.append(t.SubJet_phi.__getitem__(j))
				sub_btag.append(t.SubJet_btagDeepB.__getitem__(j))




			# =================================================================
			# FatJet Selection and SubJet Selection
			nFat = 0
			flag_Fat = 0
			dR_Fat = 0
			dR_Fat2 = 0
			dR_Sub = 0
			for j in range(t.nFatJet):

				fat_mass.append(t.FatJet_mass.__getitem__(j))
				fat_pt.append(t.FatJet_pt.__getitem__(j))
				fat_eta.append(t.FatJet_eta.__getitem__(j))
				fat_phi.append(t.FatJet_phi.__getitem__(j))
				fat_btag.append(t.FatJet_btagDeepB.__getitem__(j))
				fat_jetid.append(t.FatJet_jetId.__getitem__(j))
				fat_sdm.append(t.FatJet_msoftdrop.__getitem__(j))
				fat_tau1.append(t.FatJet_tau1.__getitem__(j))
				fat_tau2.append(t.FatJet_tau2.__getitem__(j))

				count_tFat = count_tFat + t.nFatJet

			# =================================================================
			# Resolved Jet selection
			nThin = 0
			flag_thin = 0
			dR_Res = 0
			dR_Res2 = 0

                        nThin_boost = 0
                        nThin_res = 0

			tBtag_res = 0 # added by 240502


			LED_JEC = []

			for j in range(t.nJet):

				JEC1, JEC2, JEC3, JEC4, JEC5, JEC6, JEC7, JEC8, JEC9, JEC10, JEC11, JEC12 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

                                # Jet Energy correction applied by 240405

                                if (t.Jet_pt.__getitem__(j) > 40 and t.Jet_pt.__getitem__(j) < 50):

                                        if (abs(t.Jet_eta.__getitem__(j)) < 2.5 and t.Jet_puIdDisc.__getitem__(j) > -0.42):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC1 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.5 and abs(t.Jet_eta.__getitem__(j)) < 2.75 and t.Jet_puIdDisc.__getitem__(j) > -0.09):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC2 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.75 and abs(t.Jet_eta.__getitem__(j)) < 3.0 and t.Jet_puIdDisc.__getitem__(j) > -0.14):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC3 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 3.0 and abs(t.Jet_eta.__getitem__(j)) < 5.0 and t.Jet_puIdDisc.__getitem__(j) > -0.02):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC4 = 1

					if (JEC1 == 1 or JEC2 == 1 or JEC3 == 1 or JEC4 == 1): LED_JEC.append(1)
					else: LED_JEC.append(0)

                                if (t.Jet_pt.__getitem__(j) > 30 and t.Jet_pt.__getitem__(j) < 40):
                                        if (abs(t.Jet_eta.__getitem__(j)) < 2.5 and t.Jet_puIdDisc.__getitem__(j) > -0.71):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC5 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.5 and abs(t.Jet_eta.__getitem__(j)) < 2.75 and t.Jet_puIdDisc.__getitem__(j) > -0.36):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC6 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.75 and abs(t.Jet_eta.__getitem__(j)) < 3.0 and t.Jet_puIdDisc.__getitem__(j) > -0.29):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC7 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 3.0 and abs(t.Jet_eta.__getitem__(j)) < 5.0 and t.Jet_puIdDisc.__getitem__(j) > -0.23):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC8 = 1

					if (JEC5 == 1 or JEC6 == 1 or JEC7 == 1 or JEC8 == 1): LED_JEC.append(1)
					else: LED_JEC.append(0)


                                if (t.Jet_pt.__getitem__(j) > 20 and t.Jet_pt.__getitem__(j) < 30):
                                        if (abs(t.Jet_eta.__getitem__(j)) < 2.5 and t.Jet_puIdDisc.__getitem__(j) > -0.90):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC9 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.5 and abs(t.Jet_eta.__getitem__(j)) < 2.75 and t.Jet_puIdDisc.__getitem__(j) > -0.57):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC10 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.75 and abs(t.Jet_eta.__getitem__(j)) < 3.0 and t.Jet_puIdDisc.__getitem__(j) > -0.43):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC11 = 1


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 3.0 and abs(t.Jet_eta.__getitem__(j)) < 5.0 and t.Jet_puIdDisc.__getitem__(j) > -0.42):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))
						JEC12 = 1

					if (JEC9 == 1 or JEC10 == 1 or JEC11 == 1 or JEC12 == 1): LED_JEC.append(1)
					else: LED_JEC.append(0)


                                if (t.Jet_pt.__getitem__(j) > 10 and t.Jet_pt.__getitem__(j) < 20):

					LED_JEC.append(0)
                                        if (abs(t.Jet_eta.__getitem__(j)) < 2.5 and t.Jet_puIdDisc.__getitem__(j) > -0.95):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.5 and abs(t.Jet_eta.__getitem__(j)) < 2.75 and t.Jet_puIdDisc.__getitem__(j) > -0.70):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 2.75 and abs(t.Jet_eta.__getitem__(j)) < 3.0 and t.Jet_puIdDisc.__getitem__(j) > -0.52):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))


                                        if (abs(t.Jet_eta.__getitem__(j)) >= 3.0 and abs(t.Jet_eta.__getitem__(j)) < 5.0 and t.Jet_puIdDisc.__getitem__(j) > -0.49):
                                                thin_mass.append(t.Jet_mass.__getitem__(j))
                                                thin_pt.append(t.Jet_pt.__getitem__(j))
                                                thin_eta.append(t.Jet_eta.__getitem__(j))
                                                thin_phi.append(t.Jet_phi.__getitem__(j))
                                                thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                                thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                                thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                                thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

						thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))


				if (t.Jet_pt.__getitem__(j) <= 10): LED_JEC.append(0)
				if (t.Jet_pt.__getitem__(j) == 20 or t.Jet_pt.__getitem__(j) == 30 or t.Jet_pt.__getitem__(j) == 40 or t.Jet_pt.__getitem__(j) == 50): LED_JEC.append(0)

                                if (t.Jet_pt.__getitem__(j) > 50):

					LED_JEC.append(0)
                                        thin_mass.append(t.Jet_mass.__getitem__(j))
                                        thin_pt.append(t.Jet_pt.__getitem__(j))
                                        thin_eta.append(t.Jet_eta.__getitem__(j))
                                        thin_phi.append(t.Jet_phi.__getitem__(j))
                                        thin_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
                                        thin_jetid.append(t.Jet_jetId.__getitem__(j))

                                        thin_nMu.append(t.Jet_nMuons.__getitem__(j))
                                        thin_nEl.append(t.Jet_nElectrons.__getitem__(j))

					thin_tbtag.append(t.Jet_btagDeepB.__getitem__(j))

                        # =================================================================



			# =================================================================
			cand_res_pt1, cand_res_eta1, cand_res_phi1, cand_res_mass1 = [], [], [], []
			cand_res_btag1, cand_res_jetid1, cand_res_nMu1, cand_res_nEl1 = [], [], [], []
			cand_res_tbtag1 = [] # added by 240513
			cand_res_JEC1 = [] # added by 240826

			overlap_pt, overlap_eta, overlap_phi, overlap_mass = 0, 0, 0, 0

			flag_dR_mu, flag_dR_el = 2, 2 ## meaning of numbers: 2 = default, 1 = veto, 0 = pass
			dR_res_mu, dR_res_el = 0, 0


                        btag_4jets = 0
                        non_4jets = 0

			if (hlt_mu == True and hlt_el == True): continue
			if (hlt_mu == False and hlt_el == False): continue
			if (met_filter == False): continue

			for j in range(len(thin_pt)):

				if (thin_pt[j] > 25 and abs(thin_eta[j]) < 2.4 and thin_jetid[j] >= 7):

					cand_res_pt1.append(thin_pt[j])
					cand_res_eta1.append(thin_eta[j])
					cand_res_phi1.append(thin_phi[j])
					cand_res_mass1.append(thin_mass[j])
					cand_res_btag1.append(thin_btag[j])
					cand_res_jetid1.append(thin_jetid[j])

					cand_res_JEC1.append(LED_JEC[j])
					cand_res_tbtag1.append(thin_tbtag[j])

					overlap_pt = thin_pt[j]
					overlap_eta = thin_eta[j]
					overlap_phi = thin_phi[j]
					operlap_mass = thin_mass[j]


			if (len(cand_res_pt1) < 3): continue
						
			# =================================================================


			lep1_p4 = TLorentzVector() ## for resolved Higgs

			if (flag_Mu == 1 and flag_El == 0):
				lep1_p4.SetPtEtaPhiM(llep_pt, llep_eta, llep_phi, llep_mass)
			if (flag_Mu == 0 and flag_El == 1):
				lep1_p4.SetPtEtaPhiM(el_pt, el_eta, el_phi, el_mass)


			met_p4 = TLorentzVector()
			met_p4.SetPtEtaPhiM(met_pt,0,met_phi,0)

			w_met1 = TLorentzVector()
			if (lep1_p4.Pt() > 0 and met_p4.Pt() > 0): w_met1.SetPtEtaPhiM( (lep1_p4 + met_p4).Pt(), (lep1_p4 + met_p4).Eta(), (lep1_p4 + met_p4).Phi(), (lep1_p4 + met_p4).M() )

			# =================================================================
			

			chi_h = [] ## for applying chi square method to HH reconstruction
			chi_hw = []

			jj1, jj2 = TLorentzVector(), TLorentzVector()
			bb1, bb2 = TLorentzVector(), TLorentzVector()

			b_cand_pt, b_cand_eta, b_cand_phi, b_cand_mass = [], [], [], []
			j_cand_pt, j_cand_eta, j_cand_phi, j_cand_mass = [], [], [], []

			H_chi, W_chi1, W_chi2 = [], [], []
			W_lv = TLorentzVector()

			sel_pt, sel_eta, sel_phi, sel_mass, sel_btag = [], [], [], [], []
			sel_bpt, sel_beta, sel_bphi, sel_bmass = [], [], [], []

			b_check, j_check, had_check = 0, 0, 0
			for j in range(len(cand_res_pt1)):
				if (cand_res_pt1[j] > 25 and abs(cand_res_eta1[j]) < 2.4 and cand_res_jetid1[j] >= 7):
					j_check += 1
					if (cand_res_btag1[j] > 0.3093): 
						b_check += 1
						sel_bpt.append(cand_res_pt1[j])
						sel_beta.append(cand_res_eta1[j])
						sel_bphi.append(cand_res_phi1[j])
						sel_bmass.append(cand_res_mass1[j])

			if (b_check < 2): continue

			for j in range(len(cand_res_pt1)):
				if (cand_res_pt1[j] > 25 and abs(cand_res_eta1[j]) < 2.4 and cand_res_jetid1[j] >= 7):
					if (cand_res_btag1[j] > 0.3093): continue
					had_check += 1
					sel_pt.append(cand_res_pt1[j])
					sel_eta.append(cand_res_eta1[j])
					sel_phi.append(cand_res_phi1[j])
					sel_mass.append(cand_res_mass1[j])
					sel_btag.append(cand_res_btag1[j])


			J_higgs1, J_higgs2, J_higgs3, J_higgs4 = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()
			J_higgs5 = TLorentzVector()

			BJet_1, BJet_2 = TLorentzVector(), TLorentzVector()
			Jet_1, Jet_2 = TLorentzVector(), TLorentzVector()

			BJet_1.SetPtEtaPhiM(sel_bpt[0], sel_beta[0], sel_bphi[0], sel_bmass[0])
			BJet_2.SetPtEtaPhiM(sel_bpt[1], sel_beta[1], sel_bphi[1], sel_bmass[1])

			W_higgs1 = J_higgs1 + J_higgs2
			W_higgs2 = TLorentzVector()
			Chi_had1, Chi_had2 = [], []
			Chi_W1_1, Chi_W1_2, Chi_W1_3 = [], [], []
			Chi_W2_1, Chi_W2_2, Chi_W2_3 = [], [], []

			hadronic_pt1, hadronic_eta1, hadronic_phi1, hadronic_mass1 = [], [], [], []
			hadronic_pt2, hadronic_eta2, hadronic_phi2, hadronic_mass2 = [], [], [], []
			hadronic_pt3, hadronic_eta3, hadronic_phi3, hadronic_mass3 = [], [], [], []
			hadronic_btag1, hadronic_btag2, hadronic_btag3 = [], [], []

			had_pt1, had_eta1, had_phi1, had_mass1 = [], [], [], []
			had_pt2, had_eta2, had_phi2, had_mass2 = [], [], [], []


			W_had, Had1_p4, Had2_p4 = TLorentzVector(), TLorentzVector(), TLorentzVector()
			idx_1, idx_2 = [], []
			dR_jj = []


			if (had_check == 1):
				if (abs(z_mass - sel_mass[0]) <= 10): continue
				Jet_1.SetPtEtaPhiM(sel_pt[0], sel_eta[0], sel_phi[0], sel_mass[0])
				Jet_2.SetPtEtaPhiM(0, 0, 0, 0)

			if (had_check == 2):
				Jet_1.SetPtEtaPhiM(sel_pt[0], sel_eta[0], sel_phi[0], sel_mass[0])
				Jet_2.SetPtEtaPhiM(sel_pt[1], sel_eta[1], sel_phi[1], sel_mass[1])
				if (abs(z_mass - (Jet_1 + Jet_2).M()) <= 10): continue

			if (had_check > 2):
				for j in range(len(sel_pt)):
					if (j >= len(sel_pt) - 1): continue
					J_higgs3.SetPtEtaPhiM(sel_pt[j], sel_eta[j], sel_phi[j], sel_mass[j])
					for k in range(len(sel_pt)):
						if (k <= j): continue
						J_higgs4.SetPtEtaPhiM(sel_pt[k], sel_eta[k], sel_phi[k], sel_mass[k])

						if (abs(z_mass - (J_higgs3 + J_higgs4).M()) <= 10): continue

						hadronic_pt1.append(sel_pt[j])
						hadronic_eta1.append(sel_eta[j])
						hadronic_phi1.append(sel_phi[j])
						hadronic_mass1.append(sel_mass[j])
						idx_1.append(j)

						hadronic_pt2.append(sel_pt[k])
						hadronic_eta2.append(sel_eta[k])
						hadronic_phi2.append(sel_phi[k])
						hadronic_mass2.append(sel_mass[k])
						idx_2.append(k)

						W_higgs2 = J_higgs3 + J_higgs4
						had_pt1.append(W_higgs2.Pt())
						had_eta1.append(W_higgs2.Eta())
						had_phi1.append(W_higgs2.Phi())
						had_mass1.append(W_higgs2.M())
						
						dR_jj.append(J_higgs3.DeltaR(J_higgs4))

				for j in range(len(dR_jj)):
					if (len(hadronic_pt1) == 0 or len(hadronic_pt2) == 0): continue
					if (min(dR_jj) == dR_jj[j]):
						Jet_1.SetPtEtaPhiM(hadronic_pt1[j], hadronic_eta1[j], hadronic_phi1[j], hadronic_mass1[j])
						Jet_2.SetPtEtaPhiM(hadronic_pt2[j], hadronic_eta2[j], hadronic_phi2[j], hadronic_mass2[j])
#				if (abs(z_mass - (Jet_1 + Jet_2).M()) <= 10): continue


			# ---------- Event Veto section (HH)----------

			if (met_pt == 0): continue
			if (llep_pt == 0 and el_pt == 0): continue
			if (llep_pt > 0 and el_pt > 0): continue
			if (lep1_p4.Pt() == 0): continue
			if (sel_bpt[0] == 0 or sel_bpt[1] == 0): continue
			if (Jet_1.Pt() == 0 and Jet_2.Pt() == 0): continue
#			if (W_had.Pt() == 0): continue
			if (had_check < 1): continue

			# ----------------------------------------

			#information: t_res1 = leading jet of Bjet included, t_res2 = subleading jet of Bjet included

			passed_event += 1

			# Save CSV file	
			data.append([BJet_1.Pt(),BJet_1.Eta(),BJet_1.Phi(),BJet_1.M(),
					BJet_2.Pt(),BJet_2.Eta(),BJet_2.Phi(),BJet_2.M(),
					Jet_1.Pt(),Jet_1.Eta(),Jet_1.Phi(),Jet_1.M(),
					Jet_2.Pt(),Jet_2.Eta(),Jet_2.Phi(),Jet_2.M(),
					lep1_p4.Pt(),lep1_p4.Eta(),lep1_p4.Phi(),lep1_p4.M(),
					met_p4.Pt(), met_p4.Eta(), met_p4.Phi(), met_p4.M(),
					passed_event])

        df = pd.DataFrame(data, columns= ['bjet1_mass','bjet1_pt','bjet1_eta','bjet1_phi',
					'bjet2_mass','bjet2_pt','bjet2_eta','bjet2_phi',
					'jet1_mass','jet1_pt','jet1_eta','jet1_phi',
					'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
					'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
					'MET_pT','MET_eta','MET_phi','MET_mass',
					'Event_Number'])
				

		
        df.to_csv("HH_Object_Selection.csv", header=True, index=False)

        print("Have processed %i files..."%(ii))
        print("Have processed total %i events..."%(ee))
        print("Event number finally passed is %i..."%(passed_event))

	print("")
	print("")
        print("***********************************************************")
	now2 = datetime.datetime.now()
	print("Analysis started time: ", now1.strftime("%Y-%m-%d %H:%M:%S"))
        print("Analysis finished time: ", now2.strftime("%Y-%m-%d %H:%M:%S"))
	print("")
	print("Total Pased Events: %i"%(passed_event))
        print("***********************************************************")
        print("")
        print("")



	input("Press Enter to continue...")		


if __name__ == "__main__":
	main()
