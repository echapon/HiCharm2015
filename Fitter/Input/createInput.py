#!/usr/bin/python

import os,sys
from shutil import copyfile

##########################
## INPUT PARAMETERS HERE #
##########################

bins=[
        "0.0-1.6;6.5-30.0;0.0-100.0",
        "0.0-1.6;6.5-8;0.0-100.0",
        "0.0-1.6;8-9.5;0.0-100.0",
        "0.0-1.6;9.5-30;0.0-100.0"
        ]

# provide the signal model for pp and PbPb
fcn_signal_pp="DoubleCrystalBall"
fcn_signal_pbpb="DoubleCrystalBall"

# provide the background model for pp and PbPb
fcn_bkg_pp="SecondOrderChebychev"
fcn_bkg_pbpb="SecondOrderChebychev"

# provide the name of the parameters of your models
parnames_signal_pp="N_Jpsi_PbPb;sigma1_Jpsi_PbPb;sigma2_Jpsi_PbPb;m_Jpsi_PbPb;alpha_Jpsi_PbPb;n_Jpsi_PbPb;f_Jpsi_PbPb"
parnames_signal_pbpb="N_Jpsi_PbPb;sigma1_Jpsi_PbPb;sigma2_Jpsi_PbPb;m_Jpsi_PbPb;alpha_Jpsi_PbPb;n_Jpsi_PbPb;f_Jpsi_PbPb"
parnames_bkg_pp="N_Bkg_PbPb;lambda1_Bkg_PbPb;lambda2_Bkg_PbPb;lambda3_Bkg_PbPb;lambda4_Bkg_PbPb"
parnames_bkg_pbpb="N_Bkg_PbPb;lambda1_Bkg_PbPb;lambda2_Bkg_PbPb;lambda3_Bkg_PbPb;lambda4_Bkg_PbPb"

# provide the initial parameters
parini_signal_pp="[ 0 , 250000 ];[ 0.04 , 0.01 , 0.09 ];[ 0.02 , 0.01 , 0.07 ];[ 3.096 , 2.996 , 3.196 ];[ 2.21 , 0.5 , 3000.0 ];[ 2.0 , 0.5 , 3000.0 ];[ 0.3 , 0.0 , 1.0 ]"
parini_signal_pbpb="[ 0 , 250000 ];[ 0.04 , 0.01 , 0.09 ];[ 0.02 , 0.01 , 0.07 ];[ 3.096 , 2.996 , 3.196 ];[ 2.21 , 0.5 , 3000.0 ];[ 2.0 , 0.5 , 3000.0 ];[ 0.3 , 0.0 , 1.0 ]"
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
file_sig_pp.write('rap;pt;cent;Model_JPSI_PP;' + parnames_signal_pp + '\n')
for strbin in bins:
    file_sig_pp.write(strbin + ';' + fcn_signal_pp + ';' + parini_signal_pp + '\n')

# signal, pbpb
file_sig_pbpb = open(name + '/InitialParam_MASS_JPSI_PbPb.csv', 'w')
file_sig_pbpb.write('rap;pt;cent;Model_JPSI_PbPb;' + parnames_signal_pbpb + '\n')
for strbin in bins:
    file_sig_pbpb.write(strbin + ';' + fcn_signal_pbpb + ';' + parini_signal_pbpb + '\n')

# bkg, pp
file_sig_pp = open(name + '/InitialParam_MASS_Bkg_PP.csv', 'w')
file_sig_pp.write('rap;pt;cent;Model_Bkg_PP;' + parnames_bkg_pp + '\n')
for strbin in bins:
    file_sig_pp.write(strbin + ';' + fcn_bkg_pp + ';' + parini_bkg_pp + '\n')

# bkg, pbpb
file_sig_pbpb = open(name + '/InitialParam_MASS_Bkg_PbPb.csv', 'w')
file_sig_pbpb.write('rap;pt;cent;Model_Bkg_PbPb;' + parnames_bkg_pbpb + '\n')
for strbin in bins:
    file_sig_pbpb.write(strbin + ';' + fcn_bkg_pbpb + ';' + parini_bkg_pbpb + '\n')

# finally, copy the list of input files
copyfile('InputTrees.txt',name+'/InputTrees.txt')
