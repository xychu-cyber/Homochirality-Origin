# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 17:38:46 2020

@author: Chu Xin-Yi
"""
#Psyn_LL: Probability of L-PSU synthesize L-peptide.	0.99
#Psyn_LR/D:	Probability of L-PSU synthesize R/D-peptide.	0.01
#Psyn_RL: Probability of R-PSU synthesize L-peptide.	0.99
#Psyn_RR/D:	Probability of R-PSU synthesize R/D-peptide.	0.01
#Pbind:	Probability of PSU being bound by its synthesized L-peptide.	0.1
#Prep:	Probability of PSU replication.	0.1
#Pdeg_psu:	Probability of PSU degradation.	0.1
#Pdeg_com:	Probability of PSU-pepetide complex degradation.	0.01, 0.03, 0.05, 0.07, 0.1

w = open('C:/Users/Chu Xin-Yi/Desktop/Simulation/result.txt','w')

#initial input 
psu_R = 10000000000     # R-PSU: Peptide synthesis unit don’t distinguish substrate chirality.
psu_L = 1               # L-PSU: Peptide synthesis unit prefer L-amino acids as substrates.
comx_R = 0              # R-complex: Complex of R-PSU and L-peptide.
comx_L = 0              # L-complex: Complex of L-PSU and L-peptide.
pep_R = 0               # R/D-peptide:Peptides composed of mixed L/D-amino acids or entirely D-amino acids.
pep_L = 0               # L-peptide: Peptides composed entirely of L-amino acids.

#Records of every step  
rpsu_R = []
rpsu_L = []
rpep_R = []
rpep_L = []
ee_psu = []             # (psu_L)/(psu_L+psu_R)
ee_pep = []             # (pep_L)/(pep_L+pep_R)

# simulation steps
i = 1000                # number of simulation steps 
while i > 0:
    # peptide synthesis
    pep_R += psu_R*0.99 + psu_L*0.05               # R/D-pep = R-PSU*Psyn_RR/D + L-PSU*Psyn_LR/D
    pep_L += psu_R*0.01 + psu_L*0.495              # L-pep = R-PSU*Psyn_RL + L-PSU*Psyn_LL
    # protective peptide binding
    comx_R += psu_R*0.01*0.1                       # R-complex = R-PSU*Psyn_RL*Pbind
    comx_L += psu_L*0.495*0.1                      # L-complex = L-PSU*Psyn_LL*Pbind
    # PSU degradation    
    psu_R -= (psu_R - comx_R)*0.1 + comx_R*0.01    # R-PSU = R-PSU - (R-PSU - R-complex)*Pdeg_psu -  R-complex*Pdeg_com
    psu_L -= (psu_L - comx_L)*0.1 + comx_L*0.01    # L-PSU = L-PSU - (L-PSU - L-complex)*Pdeg_psu -  L-complex*Pdeg_com
    # PSU replication
    psu_R += psu_R*0.1                             # R-PSU = R-PSU + R-PSU*Prep
    psu_L += psu_L*0.1                             # L-PSU = L-PSU + L-PSU*Prep
    # record     
    rpsu_R.append(psu_R)
    rpsu_L.append(psu_L)
    rpep_R.append(pep_R)
    rpep_L.append(pep_L)
    rpep_R.append(pep_R)
    ee_psu.append(str((psu_L)/(psu_L+psu_R)))      # the ratio of L-PSU in all PSU
    ee_pep.append(str((pep_L)/(pep_L+pep_R)))      # the ratio of L-pep in all peptide
    i -= 1
    
#print((psu_L)/(psu_L+psu_R),(pep_L)/(pep_L+pep_R))
# output
w.write('\t'.join(ee_psu)+'\n'+'\t'.join(ee_pep))
w.flush();w.close()
    