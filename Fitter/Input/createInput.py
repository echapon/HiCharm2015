#!/usr/bin/python

import os,sys
from shutil import copyfile

##########################
## INPUT PARAMETERS HERE #
##########################

bins=[
        "0.0-1.6;6.5-30.0;0.0-100.0",
        ]

# provide the signal model for pp and PbPb
fcn_signal_pp=["DoubleCrystalBall","GaussianAndCrystalBall"]
fcn_signal_pbpb=["DoubleCrystalBall","GaussianAndCrystalBall"]
fcn_signal2_pp=["DoubleCrystalBall","GaussianAndCrystalBall"]
fcn_signal2_pbpb=["DoubleCrystalBall","GaussianAndCrystalBall"]

# provide the background model for pp and PbPb
fcn_bkg_pp=["SecondOrderChebychev","SecondOrderPolynomial"]
fcn_bkg_pbpb=["SecondOrderChebychev","SecondOrderPolynomial"]

# provide the name of the parameters of your models
parnames_signal_pp="N_Jpsi_PP;sigma1_Jpsi_PP;sigma2_Jpsi_PP;m_Jpsi_PP;alpha_Jpsi_PP;n_Jpsi_PP;f_Jpsi_PP"
parnames_signal_pbpb="N_Jpsi_PbPb;sigma1_Jpsi_PbPb;sigma2_Jpsi_PbPb;m_Jpsi_PbPb;alpha_Jpsi_PbPb;n_Jpsi_PbPb;f_Jpsi_PbPb"
parnames_signal2_pp="N_Psi2S_PP;sigma1_Psi2S_PP;sigma2_Psi2S_PP;m_Psi2S_PP;alpha_Psi2S_PP;n_Psi2S_PP;f_Psi2S_PP"
parnames_signal2_pbpb="N_Psi2S_PbPb;sigma1_Psi2S_PbPb;sigma2_Psi2S_PbPb;m_Psi2S_PbPb;alpha_Psi2S_PbPb;n_Psi2S_PbPb;f_Psi2S_PbPb"
parnames_bkg_pp="N_Bkg_PP;lambda1_Bkg_PP;lambda2_Bkg_PP;lambda3_Bkg_PP;lambda4_Bkg_PP"
parnames_bkg_pbpb="N_Bkg_PbPb;lambda1_Bkg_PbPb;lambda2_Bkg_PbPb;lambda3_Bkg_PbPb;lambda4_Bkg_PbPb"

# provide the initial parameters
parini_signal_pp="[ 0 , 250000 ];[ 0.04 , 0.01 , 0.09 ];[ 0.02 , 0.01 , 0.07 ];[ 3.096 , 2.996 , 3.196 ];[ 2.21 , 0.5 , 10.0 ];[ 2.0 , 0.5 , 20. ];[ 0.3 , 0.0 , 1.0 ]"
parini_signal_pbpb="[ 0 , 250000 ];[ 0.04 , 0.01 , 0.09 ];[ 0.02 , 0.01 , 0.07 ];[ 3.096 , 2.996 , 3.196 ];[ 2.21 , 0.5 , 10.0 ];[ 2.0 , 0.5 , 20. ];[ 0.3 , 0.0 , 1.0 ]"
parini_signal2_pp=";;;;;;"
parini_signal2_pbpb=";;;;;;"
parini_bkg_pp="[0,300000];[0.05,-1.0,1.0];[0.05,-1.0,1.0];;"
parini_bkg_pbpb="[ 0 , 250000 ];[ 0.05 , -1.0 , 1.0 ];[ 0.05 , -1.0 , 1.0 ];;"

########################
# END INPUT PARAMETERS #
########################

if len(sys.argv) != 2:
   sys.exit("Please provide the name of your input (i.e. the name of the directory that I will create)")
    

name=str(sys.argv[1])

# provide the desired list of bins
# format: ymin-ymax;ptmin-ptmax;centmin-centmax

print 'I will create the input', name

os.mkdir(name)

# signal, pp
file_sig_pp = open(name + '/InitialParam_MASS_JPSI_PP.csv', 'w')
file_sig_pp.write('rap;pt;cent;Model_Jpsi_PP;' + parnames_signal_pp + '\n')
for strbin in bins:
  for strsig in fcn_signal_pp:
    file_sig_pp.write(strbin + ';' + strsig + ';' + parini_signal_pp + '\n')

# signal, pbpb
file_sig_pbpb = open(name + '/InitialParam_MASS_JPSI_PbPb.csv', 'w')
file_sig_pbpb.write('rap;pt;cent;Model_Jpsi_PbPb;' + parnames_signal_pbpb + '\n')
for strbin in bins:
  for strsig in fcn_signal_pbpb:
    file_sig_pbpb.write(strbin + ';' + strsig + ';' + parini_signal_pbpb + '\n')

# signal2, pp
file_sig_pp = open(name + '/InitialParam_MASS_PSI2S_PP.csv', 'w')
file_sig_pp.write('rap;pt;cent;Model_Psi2S_PP;' + parnames_signal2_pp + '\n')
for strbin in bins:
  for strsig in fcn_signal2_pp:
    file_sig_pp.write(strbin + ';' + strsig + ';' + parini_signal2_pp + '\n')

# signal2, pbpb
file_sig_pbpb = open(name + '/InitialParam_MASS_PSI2S_PbPb.csv', 'w')
file_sig_pbpb.write('rap;pt;cent;Model_Psi2S_PbPb;' + parnames_signal2_pbpb + '\n')
for strbin in bins:
  for strsig in fcn_signal2_pbpb:
    file_sig_pbpb.write(strbin + ';' + strsig + ';' + parini_signal2_pbpb + '\n')

# bkg, pp
file_sig_pp = open(name + '/InitialParam_MASS_BKG_PP.csv', 'w')
file_sig_pp.write('rap;pt;cent;Model_Bkg_PP;' + parnames_bkg_pp + '\n')
for strbin in bins:
  for strbkg in fcn_bkg_pp:
    file_sig_pp.write(strbin + ';' + strbkg + ';' + parini_bkg_pp + '\n')

# bkg, pbpb
file_sig_pbpb = open(name + '/InitialParam_MASS_BKG_PbPb.csv', 'w')
file_sig_pbpb.write('rap;pt;cent;Model_Bkg_PbPb;' + parnames_bkg_pbpb + '\n')
for strbin in bins:
  for strbkg in fcn_bkg_pbpb:
    file_sig_pbpb.write(strbin + ';' + strbkg + ';' + parini_bkg_pbpb + '\n')

# finally, copy the list of input files
copyfile('InputTrees.txt',name+'/InputTrees.txt')
