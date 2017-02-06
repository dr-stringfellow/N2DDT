from __future__ import division

import ROOT
from ROOT import gROOT, TFile, TTree, TChain, gPad, gDirectory, AddressOf
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
import numpy
import os

gROOT.ProcessLine(
"struct TreeStruct {\
Float_t N2DDT;\
}")

##############################################################################
## binning
nrho  = 14;
rholo = -6;
rhohi = -1;
npt   = 4;
ptlo  = 200;
pthi  = 1000;
nn2   = 500;
#n2lo  = 0;
#n2hi = 0.4;
n2lo  = 0.05;
n2hi = 0.5;


from ROOT import TreeStruct

treestruct = TreeStruct()

def main(options,args):
	
	ifile = options.ifile
	tf = ROOT.TFile(ifile);
	tt = tf.Get("events");
	nent = int(tt.GetEntries())

        Apply = True;

	if Apply: applyTransformation(tt,ifile);

######--------------------------------------------------------------------------------------------------------
######--------------------------------------------------------------------------------------------------------
def applyTransformation(tt,ifile):

	f_h2ddt = ROOT.TFile("/home/bmaier/cms/MonoHiggs/higgstagging/N2DDT/h3_n2ddt.root");
	trans_h2ddt = f_h2ddt.Get("h2ddt");

       	nent = int(tt.GetEntries())
	
	if (0!=tt.FindBranch("N2DDT")):
		tt.SetBranchStatus("N2DDT",0);


	output = ifile + "_tmp"
	ofile = ROOT.TFile(output,"RECREATE");

	otree = tt.CloneTree()
	N2DDT = array.array( 'f', [-99.0])
	o_N2DDT = otree.Branch("N2DDT" , AddressOf(treestruct,'N2DDT'), "N2DDT/F" )

	for i in range(int(tt.GetEntries())):

		tt.GetEntry(i)

		treestruct.N2DDT = -99.

		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
			sys.stdout.flush()

		jmsd_8 = tt.fj1MSD_corr
		jpt_8  = tt.fj1Pt		
		if jmsd_8 <= 0: jmsd_8 = 0.01

		rh_8 = math.log(jmsd_8*jmsd_8/jpt_8/jpt_8)
		if not math.pow(tt.fj1ECFN_1_2_10,2.00)==0.0:
			jtN2b1sd_8 = tt.fj1ECFN_2_3_10/math.pow(tt.fj1ECFN_1_2_10,2.00)
			cur_rho_index = trans_h2ddt.GetXaxis().FindBin(rh_8);
			cur_pt_index  = trans_h2ddt.GetYaxis().FindBin(jpt_8);
			if rh_8 > trans_h2ddt.GetXaxis().GetBinUpEdge( trans_h2ddt.GetXaxis().GetNbins() ): cur_rho_index = trans_h2ddt.GetXaxis().GetNbins();
			if rh_8 < trans_h2ddt.GetXaxis().GetBinLowEdge( 1 ): cur_rho_index = 1;
			if jpt_8 > trans_h2ddt.GetYaxis().GetBinUpEdge( trans_h2ddt.GetYaxis().GetNbins() ): cur_pt_index = trans_h2ddt.GetYaxis().GetNbins();
			if jpt_8 < trans_h2ddt.GetYaxis().GetBinLowEdge( 1 ): cur_pt_index = 1;
			
			treestruct.N2DDT = jtN2b1sd_8 - trans_h2ddt.GetBinContent(cur_rho_index,cur_pt_index);
		#print treestruct.N2DDT


		o_N2DDT.Fill()
		
		
	ofile.Write()
	ofile.Close()
	os.system("mv -f %s %s" % (output, ifile))


##----##----##----##----##----##----##
if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-i','--ifile', dest='ifile', default = 'file.root',help='MC/data file to add the N2DDT branch to', metavar='ifile')
	#parser.add_option('-ddt','--ddtfile', dest='ddtfile', default = '/home/bmaier/cms/MonoHiggs/higgstagging/N2DDT/h3_n2ddt.root',help='n2ddr.root file', metavar='ddtfile')

	(options, args) = parser.parse_args()
	 
	main(options,args)





