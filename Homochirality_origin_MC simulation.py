# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 17:38:46 2020
@author: Chu Xin-Yi
"""
w = open('result.txt','w')

#initial input 
psu_R = 10**10          # R-PSU: Peptide synthesis unit doesn’t distinguish substrate chirality.
psu_D = 1               # D-PSU: Peptide synthesis unit prefers D-amino acids as substrates.
psu_L = 1               # L-PSU: Peptide synthesis unit prefers L-amino acids as substrates.
comx_R_R = 0            # Complex of R-PSU and R-peptide.
comx_R_D = 0            # Complex of R-PSU and D-peptide.
comx_R_L = 0            # Complex of R-PSU and L-peptide.
comx_D_D = 0            # Complex of D-PSU and D-peptide.
comx_L_L = 0            # Complex of L-PSU and L-peptide.
pep_R = 0               # R-peptide:Peptides composed of mixed L/D-amino acids.
pep_D = 0               # D-peptide: Peptides composed entirely of L-amino acids.
pep_L = 0               # L-peptide: Peptides composed entirely of L-amino acids.

Psyn_RR = 0.99	         # Probability of R-PSU synthesize R-peptide.
Psyn_RD = 0.005	      # Probability of R-PSU synthesize D-peptide. Psyn_RD = (1-Psyn_RR)*0.5
Psyn_RL = 0.005	      # Probability of R-PSU synthesize L-peptide. Psyn_RL = (1-Psyn_RR)*0.5
Psyn_DD = 0.5	         # Probability of D-PSU synthesize D-peptide.
Psyn_LL = 0.5	         # Probability of L-PSU synthesize L-peptide.
Pbind_R = 0.01	         # Probability of PSU being bound by its synthesized R-peptide.
Pbind_D = 0.01	         # Probability of PSU being bound by its synthesized D-peptide.
Pbind_L = 0.03          # Probability of PSU being bound by its synthesized L-peptide.
Pdeg_psu = 0.5	         # Probability of PSU degradation.
Pdeg_com = 0.25	      # Probability of PSU-pepetide complex degradation.
Prep = 1	               # Probability of PSU replication.

#Records of every step  
rpsu_R = []
rpsu_D = []
rpsu_L = []
rpep_R = []
rpep_D = []
rpep_L = []
ee_psu = []             # ratio of L-PSU: (psu_L)/(psu_L+psu_R+psu_D)
ee_pep = []             # ratio of L-peptide: (pep_L)/(pep_L+pep_R+psu_R)

# simulation steps
i = 12000                 # number of simulation steps 
while i > 0:
    # peptide synthesis
    pep_R += psu_R * Psyn_RR                       
    pep_D += psu_R * Psyn_RD + psu_D * Psyn_DD        
    pep_L += psu_R * Psyn_RL + psu_L * Psyn_LL        
    
    # protective peptide binding
    comx_R_R = psu_R * Psyn_RR * Pbind_R          
    comx_R_D = psu_R * Psyn_RD * Pbind_D            
    comx_R_L = psu_R * Psyn_RL * Pbind_L            
    comx_D_D = psu_D * Psyn_DD * Pbind_D                     
    comx_L_L = psu_L * Psyn_LL * Pbind_L                     
    
    # PSU degradation    
    psu_R -= (psu_R - comx_R_R - comx_R_D - comx_R_L) * Pdeg_psu + (comx_R_R + comx_R_D + comx_R_L) * Pdeg_com    
    psu_D -= (psu_D - comx_D_D) * Pdeg_psu + comx_D_D * Pdeg_com    
    psu_L -= (psu_L - comx_L_L) * Pdeg_psu + comx_L_L * Pdeg_com    
    
    # PSU replication
    psu_R += psu_R * Prep                           
    psu_D += psu_D * Prep                           
    psu_L += psu_L * Prep                           
    
    # record     
    rpsu_R.append(str(psu_R))
    rpsu_D.append(str(psu_D))
    rpsu_L.append(str(psu_L))
    rpep_R.append(str(pep_R))
    rpep_D.append(str(pep_D))
    rpep_L.append(str(pep_L))
    ee_psu.append(str((psu_L)/(psu_L + psu_R + psu_D)))      
    ee_pep.append(str((pep_L)/(pep_L + pep_R + pep_D)))      
    i -= 1
    
#print((psu_L)/(psu_L+psu_R),(pep_L)/(pep_L+pep_R))
# output
w.write('\t'.join(rpsu_R)+'\n'+'\t'.join(rpsu_D)+'\n'+'\t'.join(rpsu_L)+'\n'+'\t'.join(rpep_R)+'\n'+'\t'.join(rpep_D)+'\n'+'\t'.join(rpep_L)+'\n'+'\t'.join(ee_psu)+'\n'+'\t'.join(ee_pep))
w.flush();w.close()
