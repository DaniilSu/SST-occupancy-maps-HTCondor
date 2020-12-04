#!/usr/bin/env python2
import argparse
import numpy as np
import ROOT as r
# Fix https://root-forum.cern.ch/t/pyroot-hijacks-help/15207 :
r.PyConfig.IgnoreCommandLineOptions = True
import shipunit as u
import rootUtils as ut
from decorators import *
import shipRoot_conf
shipRoot_conf.configure()

#inputdir = "/eos/experiment/ship/user/dasukhon/"
r.gErrorIgnoreLevel = r.kWarning
r.gROOT.SetBatch(True)
parser = argparse.ArgumentParser(description='Script to create flux maps.')
parser.add_argument(
	'inputfile',
	help='''Simulation results to use as input. '''
	'''Supports retrieving files from EOS via the XRootD protocol.''')
parser.add_argument(
	'geofile',
	help='''Geometry file to use. '''
	'''Supports retrieving files from EOS via the XRootD protocol.''')
parser.add_argument(
	'basedir',
	help='''Base directory which contains simulation and geometry files. '''
	'''Default location /eos/experiment/ship/user/dasukhon/ ''')
parser.add_argument(
	'-o',
	'--outputfile',
	default='magnet_maps.root',
	help='''File to write the flux maps to. '''
	'''Will be recreated if it already exists.''')
parser.add_argument(
	'-r',
	'--ratesfile',
	default='hit_rates.txt',
	help='''File to write the hit rates per cm^2 * t_spill. '''
	'''Will be recreated if it already exists.''')
args = parser.parse_args()
g = r.TFile.Open(args.geofile, 'read')
sGeo = g.FAIRGeom
f = r.TFile.Open(args.outputfile, 'recreate')
f.cd()
maxpt = 10. * u.GeV
maxp = 360. * u.GeV
h = {}

# Define histograms
ut.bookHist(h, 'magyoke_pid', 'Particle ID;a.u.', 1200, -600, 600)
ut.bookHist(h, 'CV_pid', 'Particle ID;a.u.', 1200, -600, 600)
for station in range(1, 5):
	ut.bookHist(h, 'T{}_pid'.format(station), 'Particle ID;a.u.', 1200, -600, 600)
	ut.bookHist(h, 'MCoil{}_pid'.format(station), 'Particle ID;a.u.', 1200, -600, 600)
for suffix, title in [('pos_mu', '#mu^{+} hits'), ('neg_mu', '#mu^{-} hits'), ('all', 'All hits')]:
	if 'mu' in suffix:
		for prefix in ['entry', 'exit']:
			ut.bookHist(h, 'CV_{}_{}'.format(prefix,suffix),
							'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32, -720,
							720)
			ut.bookHist(h, 'magyoke_{}_{}'.format(prefix,suffix),
							'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32, -720,
							720)
			ut.bookHist(h, 'SHiPMagnet_{}_{}'.format(prefix,suffix),
							'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32, -720,
							720)
		for prefix in ['entry', 'exit']:
			ut.bookHist(h, 'MCoil_{}_{}'.format(prefix,suffix),
						'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32,
						-720, 720)
		for station in range(1, 5):
			ut.bookHist(h, 'T{}_{}'.format(station,suffix),
							'{};x[cm];y[cm]'.format(title), 12, -270, +270, 24,
							-540, 540)
			ut.bookHist(h, 'Time_T{}_{}'.format(station,suffix),'{};t[ns]'.format(title),150,0,1500)
			for view in ['x1','u','v','x2']:
				for plane in range (0, 2):
					ut.bookHist(h, 'T{}_{}_plane{}_{}'.format(station, view, plane, suffix), '{};x[cm];y[cm]'.format(title), 12, -270, +270, 24, -540, 540)
	else:
		for prefix in ['entry', 'exit']:
			ut.bookHist(h, 'CV_{}_{}'.format(prefix,suffix),
						'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32, -720,
						720)
			ut.bookHist(h, 'magyoke_{}_{}'.format(prefix,suffix),
						'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32, -720,
						720)
			ut.bookHist(h, 'SHiPMagnet_{}_{}'.format(prefix,suffix),
						'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32, -720,
						720)
		for prefix in ['entry', 'exit']:
			ut.bookHist(h, 'MCoil_{}_{}'.format(prefix,suffix),
					'{};x[cm];y[cm]'.format(title), 20, -450, +450, 32,
					-720, 720)
		for station in range(1, 5):
			ut.bookHist(h, 'T{}_{}'.format(station, suffix),
						'{};x[cm];y[cm]'.format(title), 12, -270, +270, 24,
						-540, 540)
			ut.bookHist(h, 'Time_T{}_{}'.format(station,suffix),'{};t[ns]'.format(title),150,0,1500)
			for view in ['x1','u','v','x2']:
				for plane in range (0, 2):
					ut.bookHist(h, 'T{}_{}_plane{}_{}'.format(station, view, plane, suffix), '{};x[cm];y[cm]'.format(title), 12, -270, +270, 24, -540, 540)
for charge, title in [('pos','#mu^{+}'),('neg','#mu^{-}')]:
	for suffix in ['', '_original']:
		ut.bookHist(h, '{}_mu_p{}'.format(charge,suffix), '{};p[GeV];'.format(title), 100, 0, maxp)
		ut.bookHist(h, '{}_mu_pt{}'.format(charge,suffix), '{};p_t[GeV];'.format(title), 100, 0,
					maxpt)
		ut.bookHist(h, '{}_mu_ppt{}'.format(charge,suffix), '{};p[GeV];p_t[GeV];'.format(title),
					100, 0, maxp, 100, 0, maxpt)
	for station in range(1, 5):
		ut.bookHist(h, '{}_mu_p_T{}_x'.format(charge,station), '{}_T{};x[cm];p[GeV];'.format(station,title),
				12, -270, +270, 20, 0, maxp)
		ut.bookHist(h, '{}_mu_p_T{}_y'.format(charge,station), '{}_T{};y[cm];p[GeV];'.format(station,title),
				24, -540, +540, 20, 0, maxp)
		ut.bookHist(h, '{}_mu_pt_T{}_x'.format(charge,station), '{}_T{};x[cm];p_t[GeV];'.format(station,title),
				12, -270, +270, 20, 0, maxpt)
		ut.bookHist(h, '{}_mu_pt_T{}_y'.format(charge,station), '{}_T{};y[cm];p_t[GeV];'.format(station,title),
				24, -540, +540, 20, 0, maxpt)
	for prefix in ['entry', 'exit']:
		ut.bookHist(h, 'CV_{}_{}_mu_p_x'.format(prefix,charge), '{};x[cm];p[GeV];'.format(title),
				20, -450, +450, 20, 0, maxp)
		ut.bookHist(h, 'CV_{}_{}_mu_p_y'.format(prefix,charge), '{};y[cm];p[GeV];'.format(title),
				32, -720, +720, 20, 0, maxp)
		ut.bookHist(h, 'CV_{}_{}_mu_pt_x'.format(prefix,charge), '{};x[cm];p_t[GeV];'.format(title),
				20, -450, +450, 20, 0, maxpt)
		ut.bookHist(h, 'CV_{}_{}_mu_pt_y'.format(prefix,charge), '{};y[cm];p_t[GeV];'.format(title),
				32, -720, +720, 20, 0, maxpt)
		ut.bookHist(h, 'magyoke_{}_{}_mu_p_x'.format(prefix,charge), '{};x[cm];p[GeV];'.format(title),
				20, -450, +450, 20, 0, maxp)
		ut.bookHist(h, 'magyoke_{}_{}_mu_p_y'.format(prefix,charge), '{};y[cm];p[GeV];'.format(title),
				32, -720, +720, 20, 0, maxp)
		ut.bookHist(h, 'magyoke_{}_{}_mu_pt_x'.format(prefix,charge), '{};x[cm];p_t[GeV];'.format(title),
				20, -450, +450, 20, 0, maxpt)
		ut.bookHist(h, 'magyoke_{}_{}_mu_pt_y'.format(prefix,charge), '{};y[cm];p_t[GeV];'.format(title),
				32, -720, +720, 20, 0, maxpt)
		ut.bookHist(h, 'MCoil_{}_{}_mu_p_x'.format(prefix,charge), '{};x[cm];p[GeV];'.format(title),
				20, -450, +450, 20, 0, maxp)
		ut.bookHist(h, 'MCoil_{}_{}_mu_p_y'.format(prefix,charge), '{};y[cm];p[GeV];'.format(title),
				32, -720, +720, 20, 0, maxp)
		ut.bookHist(h, 'MCoil_{}_{}_mu_pt_x'.format(prefix,charge), '{};x[cm];p_t[GeV];'.format(title),
				20, -450, +450, 20, 0, maxpt)
		ut.bookHist(h, 'MCoil_{}_{}_mu_pt_y'.format(prefix,charge), '{};y[cm];p_t[GeV];'.format(title),
				32, -720, +720, 20, 0, maxpt)

def viewSelect(x):
	return {
		0:'x1',
		1:'u',
		2: 'v',
		3:'x2',
	}.get(x,'x1')

ch = r.TChain("cbmsim")
ch.Add(args.inputfile)
print "Directory ", args.basedir
n = ch.GetEntries()
print "Number of entries ", n
i = 0
for event in ch:
	if i % 10000 == 0:
		print '{}/{}'.format(i, n)
	
	muon = False
	muonid = 1

	i = i + 1
	for hit in event.strawtubesPoint:
		if hit:
			if not hit.GetEnergyLoss() > 0:
				continue
			x = hit.GetX()
			y = hit.GetY()
			z = hit.GetZ()
			px = hit.GetPx()
			py = hit.GetPy()
			pz = hit.GetPz()
			pt = np.hypot(px, py)
			P = np.hypot(pz, pt)
			pid = hit.PdgCode()
			tid = hit.GetTrackID()
			if P > 0:
				assert tid > 0
				weight = event.MCTrack[tid].GetWeight()
				assert pid not in [12, -12, 14, -14, 16, -16]
				t = hit.GetTime()
				detector_ID = hit.GetDetectorID()
				station = detector_ID / 10000000
				viewnb = (detector_ID - station * 10000000) / 1000000
				view = viewSelect(viewnb)
				plane = (detector_ID - station * 10000000 - viewnb * 1000000) / 100000
				h['T{}_all'.format(station)].Fill(x, y, weight)
				h['T{}_{}_plane{}_all'.format(station, view, plane)].Fill(x, y, weight)
				h['T{}_pid'.format(station)].Fill(pid, weight)
				h['Time_T{}_all'.format(station)].Fill(t,weight)
				if abs(pid) == 13:
					muon = True
					muonid = tid
					if pid == 13:
						h['T{}_neg_mu'.format(station)].Fill(x, y, weight)
						h['T{}_{}_plane{}_neg_mu'.format(station, view, plane)].Fill(x, y, weight)
						h['Time_T{}_neg_mu'.format(station)].Fill(t,weight)
						h['neg_mu_p'].Fill(P, weight)
						h['neg_mu_pt'].Fill(pt, weight)
						h['neg_mu_p_T{}_x'.format(station)].Fill(x, P, weight)
						h['neg_mu_p_T{}_y'.format(station)].Fill(y, P, weight)
						h['neg_mu_pt_T{}_x'.format(station)].Fill(x, pt, weight)
						h['neg_mu_pt_T{}_y'.format(station)].Fill(y, pt, weight)
						h['neg_mu_ppt'].Fill(P, pt, weight)
					else:
						h['T{}_pos_mu'.format(station)].Fill(x, y, weight)
						h['T{}_{}_plane{}_pos_mu'.format(station, view, plane)].Fill(x, y, weight)
						h['Time_T{}_pos_mu'.format(station)].Fill(t,weight)
						h['pos_mu_p'].Fill(P, weight)
						h['pos_mu_pt'].Fill(pt, weight)
						h['pos_mu_p_T{}_x'.format(station)].Fill(x, P, weight)
						h['pos_mu_p_T{}_y'.format(station)].Fill(y, P, weight)
						h['pos_mu_pt_T{}_x'.format(station)].Fill(x, pt, weight)
						h['pos_mu_pt_T{}_y'.format(station)].Fill(y, pt, weight)
						h['pos_mu_ppt'].Fill(P, pt, weight)
	for hit in event.vetoPoint:
		if hit:
			if not hit.GetEnergyLoss() > 0:
				continue
			LastX = hit.LastPoint().x()
			LastY = hit.LastPoint().y()
			LastZ = hit.LastPoint().z()
			MidX = hit.GetX()
			MidY = hit.GetY()
			MidZ = hit.GetZ()
			FirstX = 2 * MidX - LastX
			FirstY = 2 * MidY - LastY
			FirstZ = 2 * MidZ - LastZ
			px = hit.GetPx()
			py = hit.GetPy()
			pz = hit.GetPz()
			pt = np.hypot(px, py)
			P = np.hypot(pz, pt)
			pid = hit.PdgCode()
			tid = hit.GetTrackID()
			detector_ID = hit.GetDetectorID()
			assert tid > 0
			weight = event.MCTrack[tid].GetWeight()
			assert pid not in [12, -12, 14, -14, 16, -16]
			detectorEntry = sGeo.FindNode(FirstX, FirstY, FirstZ).GetName()
			detectorExit = sGeo.FindNode(LastX, LastY, LastZ).GetName()
			detector = sGeo.FindNode(MidX, MidY, MidZ).GetName()
			# SHiPMagnet:
			if 'SHiPMagnet' in detector:
				h['SHiPMagnet_entry_all'].Fill(FirstX, FirstY, weight)
				h['SHiPMagnet_exit_all'].Fill(LastX, LastY, weight)
				h['SHiPMagnet_pid'].Fill(pid, weight)
				if abs(pid) == 13:
					muon = True
					muonid = tid
					if pid == 13:
						h['neg_mu_p'].Fill(P, weight)
						h['neg_mu_pt'].Fill(pt, weight)
						h['neg_mu_ppt'].Fill(P, pt, weight)
						h['SHiPMagnet_entry_neg_mu'].Fill(FirstX, FirstY, weight)
						h['SHiPMagnet_exit_neg_mu'].Fill(LastX, LastY, weight)
					else:
						h['pos_mu_p'].Fill(P, weight)
						h['pos_mu_pt'].Fill(pt, weight)
						h['pos_mu_ppt'].Fill(P, pt, weight)
						h['SHiPMagnet_entry_pos_mu'].Fill(FirstX, FirstY, weight)
						h['SHiPMagnet_exit_pos_mu'].Fill(LastX, LastY, weight)
				continue
			# Yoke of SHiPMagnet
			if 'magyoke' in detector:
				h['magyoke_entry_all'].Fill(FirstX, FirstY, weight)
				h['magyoke_exit_all'].Fill(LastX, LastY, weight)
				h['magyoke_pid'].Fill(pid, weight)
				if abs(pid) == 13:
					muon = True
					muonid = tid
					if pid == 13:
						h['neg_mu_p'].Fill(P, weight)
						h['neg_mu_pt'].Fill(pt, weight)
						h['neg_mu_ppt'].Fill(P, pt, weight)
						h['magyoke_entry_neg_mu'].Fill(FirstX, FirstY, weight)
						h['magyoke_exit_neg_mu'].Fill(LastX, LastY, weight)
						h['magyoke_entry_neg_mu_p_x'].Fill(FirstX, P, weight)
						h['magyoke_exit_neg_mu_p_x'].Fill(LastX, P, weight)
						h['magyoke_entry_neg_mu_pt_x'].Fill(FirstX, pt, weight)
						h['magyoke_exit_neg_mu_pt_x'].Fill(LastX, pt, weight)
						h['magyoke_entry_neg_mu_p_y'].Fill(FirstY, P, weight)
						h['magyoke_exit_neg_mu_p_y'].Fill(LastY, P, weight)
						h['magyoke_entry_neg_mu_pt_y'].Fill(FirstY, pt, weight)
						h['magyoke_exit_neg_mu_pt_y'].Fill(LastY, pt, weight)
					else:
						h['pos_mu_p'].Fill(P, weight)
						h['pos_mu_pt'].Fill(pt, weight)
						h['pos_mu_ppt'].Fill(P, pt, weight)
						h['magyoke_entry_pos_mu'].Fill(FirstX, FirstY, weight)
						h['magyoke_exit_pos_mu'].Fill(LastX, LastY, weight)
						h['magyoke_entry_pos_mu_p_x'].Fill(FirstX, P, weight)
						h['magyoke_exit_pos_mu_p_x'].Fill(LastX, P, weight)
						h['magyoke_entry_pos_mu_pt_x'].Fill(FirstX, pt, weight)
						h['magyoke_exit_pos_mu_pt_x'].Fill(LastX, pt, weight)
						h['magyoke_entry_pos_mu_p_y'].Fill(FirstY, P, weight)
						h['magyoke_exit_pos_mu_p_y'].Fill(LastY, P, weight)
						h['magyoke_entry_pos_mu_pt_y'].Fill(FirstY, pt, weight)
						h['magyoke_exit_pos_mu_pt_y'].Fill(LastY, pt, weight)
				continue
			if 'MCoil' in detector:
				h['MCoil_entry_all'].Fill(FirstX, FirstY, weight)
				h['MCoil_exit_all'].Fill(LastX, LastY, weight)
				if 'MCoil1' in detector:
					h['MCoil1_pid'].Fill(pid, weight)
				if 'MCoil2' in detector:
					h['MCoil2_pid'].Fill(pid, weight)
				if 'MCoil3' in detector:
					h['MCoil3_pid'].Fill(pid, weight)
				if 'MCoil4' in detector:
					h['MCoil4_pid'].Fill(pid, weight)
				if abs(pid) == 13:
					muon = True
					muonid = tid
					if pid == 13:
						h['neg_mu_p'].Fill(P, weight)
						h['neg_mu_pt'].Fill(pt, weight)
						h['neg_mu_ppt'].Fill(P, pt, weight)
						h['MCoil_entry_neg_mu'].Fill(FirstX, FirstY, weight)
						h['MCoil_exit_neg_mu'].Fill(LastX, LastY, weight)
						h['MCoil_entry_neg_mu_p_x'].Fill(FirstX, P, weight)
						h['MCoil_exit_neg_mu_p_x'].Fill(LastX, P, weight)
						h['MCoil_entry_neg_mu_pt_x'].Fill(FirstX, pt, weight)
						h['MCoil_exit_neg_mu_pt_x'].Fill(LastX, pt, weight)
						h['MCoil_entry_neg_mu_p_y'].Fill(FirstY, P, weight)
						h['MCoil_exit_neg_mu_p_y'].Fill(LastY, P, weight)
						h['MCoil_entry_neg_mu_pt_y'].Fill(FirstY, pt, weight)
						h['MCoil_exit_neg_mu_pt_y'].Fill(LastY, pt, weight)
					else:
						h['pos_mu_p'].Fill(P, weight)
						h['pos_mu_pt'].Fill(pt, weight)
						h['pos_mu_ppt'].Fill(P, pt, weight)
						h['MCoil_entry_pos_mu'].Fill(FirstX, FirstY, weight)
						h['MCoil_exit_pos_mu'].Fill(LastX, LastY, weight)
						h['MCoil_entry_pos_mu_p_x'].Fill(FirstX, P, weight)
						h['MCoil_exit_pos_mu_p_x'].Fill(LastX, P, weight)
						h['MCoil_entry_pos_mu_pt_x'].Fill(FirstX, pt, weight)
						h['MCoil_exit_pos_mu_pt_x'].Fill(LastX, pt, weight)
						h['MCoil_entry_pos_mu_p_y'].Fill(FirstY, P, weight)
						h['MCoil_exit_pos_mu_p_y'].Fill(LastY, P, weight)
						h['MCoil_entry_pos_mu_pt_y'].Fill(FirstY, pt, weight)
						h['MCoil_exit_pos_mu_pt_y'].Fill(LastY, pt, weight)
				continue
			if 'CV' in detector:
				h['CV_entry_all'].Fill(FirstX, FirstY, weight)
				h['CV_exit_all'].Fill(LastX, LastY, weight)
				h['CV_pid'].Fill(pid, weight)
				if abs(pid) == 13:
					muon = True
					muonid = tid
					if pid == 13:
						h['neg_mu_p'].Fill(P, weight)
						h['neg_mu_pt'].Fill(pt, weight)
						h['neg_mu_ppt'].Fill(P, pt, weight)
						h['CV_entry_neg_mu'].Fill(FirstX, FirstY, weight)
						h['CV_exit_neg_mu'].Fill(LastX, LastY, weight)
						h['CV_entry_neg_mu_p_x'].Fill(FirstX, P, weight)
						h['CV_exit_neg_mu_p_x'].Fill(LastX, P, weight)
						h['CV_entry_neg_mu_pt_x'].Fill(FirstX, pt, weight)
						h['CV_exit_neg_mu_pt_x'].Fill(LastX, pt, weight)
						h['CV_entry_neg_mu_p_y'].Fill(FirstY, P, weight)
						h['CV_exit_neg_mu_p_y'].Fill(LastY, P, weight)
						h['CV_entry_neg_mu_pt_y'].Fill(FirstY, pt, weight)
						h['CV_exit_neg_mu_pt_y'].Fill(LastY, pt, weight)
					else:
						h['pos_mu_p'].Fill(P, weight)
						h['pos_mu_pt'].Fill(pt, weight)
						h['pos_mu_ppt'].Fill(P, pt, weight)
						h['CV_entry_pos_mu'].Fill(FirstX, FirstY, weight)
						h['CV_exit_pos_mu'].Fill(LastX, LastY, weight)
						h['CV_entry_pos_mu_p_x'].Fill(FirstX, P, weight)
						h['CV_exit_pos_mu_p_x'].Fill(LastX, P, weight)
						h['CV_entry_pos_mu_pt_x'].Fill(FirstX, pt, weight)
						h['CV_exit_pos_mu_pt_x'].Fill(LastX, pt, weight)
						h['CV_entry_pos_mu_p_y'].Fill(FirstY, P, weight)
						h['CV_exit_pos_mu_p_y'].Fill(LastY, P, weight)
						h['CV_entry_pos_mu_pt_y'].Fill(FirstY, pt, weight)
						h['CV_exit_pos_mu_pt_y'].Fill(LastY, pt, weight)
				continue
#			print 'Unidentified vetoPoint.'
	if muon:
		weight = event.MCTrack[muonid].GetWeight()
		pid = event.MCTrack[muonid].GetPdgCode()
		if pid == 13:
			h['neg_mu_p_original'].Fill(event.MCTrack[muonid].GetP(), weight)
			h['neg_mu_pt_original'].Fill(event.MCTrack[muonid].GetPt(), weight)
			h['neg_mu_ppt_original'].Fill(event.MCTrack[muonid].GetP(),
									  event.MCTrack[muonid].GetPt(), weight)
		if pid == -13:
			h['pos_mu_p_original'].Fill(event.MCTrack[muonid].GetP(), weight)
			h['pos_mu_pt_original'].Fill(event.MCTrack[muonid].GetPt(), weight)
			h['pos_mu_ppt_original'].Fill(event.MCTrack[muonid].GetP(),
									  event.MCTrack[muonid].GetPt(), weight)
print 'Event loop done'

rates = open(args.ratesfile,'w')
for key in h:
	classname = h[key].Class().GetName()
	if 'TH' in classname or 'TP' in classname:
		if not 'pid' in key:
			h[key].SetOption("colz")
			rate = (h[key].GetMaximum()) / (45. * 45.)
			line = str(key) + ": " + str(round(rate,1)) + " hits/(cm^2*spill)\n"
			rates.write(line)
		h[key].Write()
f.Close()
rates.close()

