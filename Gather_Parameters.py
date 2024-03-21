# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:55:42 2019

@author: Sageanne Senneff
"""

def Load_Parameters(self, no_mu, filepath, i):
    
    ## Preload Motoneuron Parameters
    SomaDict = {} # Soma
    DendDict = {} # Dendrites (D1-D4)
    GenDict  = {} # General         
    
    myParams = []
    with open(filepath + '/FivecompthTA100/FivecompthTA100_%d.hoc' % i, 'rt') as myFile:
        for myParam in myFile:
            myParams.append(myParam.rstrip('\n'))
    SomaDict = {'soma.diam':float(myParams[0].split('=')[1]),'soma.L':float(myParams[1].split('=')[1]),'soma.g_pas':float(myParams[2].split('=')[1]),
            'soma.e_pas':float(myParams[3].split('=')[1]),'soma.gbar_na3rp':float(myParams[4].split('=')[1]),'soma.gbar_naps':float(myParams[5].split('=')[1]), 
            'soma.sh_na3rp':float(myParams[6].split('=')[1]),'soma.sh_naps':float(myParams[7].split('=')[1]),'soma.ar_na3rp':float(myParams[8].split('=')[1]),
            'soma.ar_naps':float(myParams[9].split('=')[1]),'soma.gMax_kdrRL':float(myParams[10].split('=')[1]),'soma.gcamax_mAHP':float(myParams[11].split('=')[1]),
            'soma.gkcamax_mAHP':float(myParams[12].split('=')[1]),'soma.taur_mAHP':float(myParams[13].split('=')[1]),'soma.ek':float(myParams[14].split('=')[1]),
            'soma.Ra':float(myParams[15].split('=')[1]),'soma.cm':float(myParams[16].split('=')[1]),'soma.ghbar_gh':float(myParams[17].split('=')[1]),'soma.half_gh':float(myParams[18].split('=')[1])}
    DendDict = {'dend.L':float(myParams[20].split('=')[1]),'dend.diam':float(myParams[21].split('=')[1]),'dend.g_pas':float(myParams[22].split('=')[1]), 
            'dend.e_pas':float(myParams[23].split('=')[1]),'dend.gcabar_L_Ca_inact':float(myParams[24].split('=')[1]),'dend.Ra':float(myParams[25].split('=')[1]), 
            'dend.cm':float(myParams[26].split('=')[1]),'dend.ghbar_gh':float(myParams[27].split('=')[1]),'dend.half_gh':float(myParams[28].split('=')[1]), 
            'd1.gcabar_L_Ca_inact':float(myParams[29+1].split('=')[1]),'d2.gcabar_L_Ca_inact':float(myParams[30+1].split('=')[1]),
            'd3.gcabar_L_Ca_inact':float(myParams[31+1].split('=')[1]),'d4.gcabar_L_Ca_inact':float(myParams[32+1].split('=')[1])}  
    GenDict = {'qinf_na3rp':float(myParams[34].split('=')[1]), 'thinf_na3rp':float(myParams[35].split('=')[1]),'vslope_naps':float(myParams[36].split('=')[1]),
            'asvh_naps':float(myParams[37].split('=')[1]),'bsvh_naps':float(myParams[38].split('=')[1]),'mvhalfca_mAHP':float(myParams[39].split('=')[1]),
            'mtauca_mAHP':float(myParams[40].split('=')[1]),'celsius':float(myParams[41].split('=')[1]),'theta_m_L_Ca_inact':float(myParams[42].split('=')[1]),
            'tau_m_L_Ca_inact':float(myParams[43].split('=')[1]),'theta_h_L_Ca_inact':float(myParams[44].split('=')[1]),'tau_h_L_Ca_inact':float(myParams[45].split('=')[1]),
            'kappa_h_L_Ca_inact':float(myParams[46].split('=')[1]),'htau_gh':float(myParams[47].split('=')[1]),'mVh_kdrRL':float(myParams[48].split('=')[1]),
            'tmin_kdrRL':float(myParams[49].split('=')[1]),'taumax_kdrRL':float(myParams[50].split('=')[1])}      
        
    return SomaDict, DendDict, GenDict

def Set_Parameters_Soma(self, no_mu, filepath, CaPIC, i):
    
    SomaDict, DendDict, GenDict = Load_Parameters(self, no_mu, filepath, i)
    
    ## Set Biophysical Parameters            
    diam = SomaDict['soma.diam'] 
    L = SomaDict['soma.L']
    Ra = SomaDict['soma.Ra']
    cm = SomaDict['soma.cm']
    ## Set Ion Conductances
    gbar_na3rp = SomaDict['soma.gbar_na3rp']
    gbar_naps = SomaDict['soma.gbar_naps']
    gMax_kdrRL = SomaDict['soma.gMax_kdrRL']
    gcamax_mAHP = SomaDict['soma.gcamax_mAHP']
    gkcamax_mAHP = SomaDict['soma.gkcamax_mAHP'] 
    ghbar_gh = SomaDict['soma.ghbar_gh']
    g_pas = SomaDict['soma.g_pas']
    ek = SomaDict['soma.ek']
    e_pas = SomaDict['soma.e_pas'] 
    ## Set Time Constants            
    taur_mAHP = SomaDict['soma.taur_mAHP']
    mtauca_mAHP = GenDict['mtauca_mAHP']
    htau_gh = GenDict['htau_gh']
    tmin_kdrRL = GenDict['tmin_kdrRL']
    taumax_kdrRL = GenDict['taumax_kdrRL']
    ## Set Voltage Parameters            
    sh_na3rp = SomaDict['soma.sh_na3rp']
    sh_naps = SomaDict['soma.sh_naps']
    ar_na3rp = SomaDict['soma.ar_na3rp']
    ar_naps = SomaDict['soma.ar_naps']
    half_gh = SomaDict['soma.half_gh'] 
    qinf_na3rp  = GenDict['qinf_na3rp']
    thinf_na3rp = GenDict['thinf_na3rp']
    vslope_naps = GenDict['vslope_naps']
    asvh_naps = GenDict['asvh_naps']
    bsvh_naps = GenDict['bsvh_naps']
    mvhalfca_mAHP = GenDict['mvhalfca_mAHP']
    mVh_kdrRL = GenDict['mVh_kdrRL']
    
    return diam, L, Ra, cm, gbar_na3rp, gbar_naps, gMax_kdrRL, gcamax_mAHP, gkcamax_mAHP, ghbar_gh, g_pas, ek, e_pas, taur_mAHP, \
mtauca_mAHP, htau_gh, tmin_kdrRL, taumax_kdrRL, sh_na3rp, sh_naps, ar_na3rp, ar_naps, half_gh, qinf_na3rp, thinf_na3rp, vslope_naps, \
asvh_naps, bsvh_naps, mvhalfca_mAHP, mVh_kdrRL
        
def Set_Parameters_D1(self, no_mu, filepath, CaPIC, i):
    
    SomaDict, DendDict, GenDict = Load_Parameters(self, no_mu, filepath, i)
        
    ## Set Biophysical Parameters   
    diam = DendDict['dend.diam']  
    L = DendDict['dend.L']
    Ra = DendDict['dend.Ra']
    cm = DendDict['dend.cm']
    gcabar_L_Ca_inact = CaPIC*(DendDict['d1.gcabar_L_Ca_inact'])
    ghbar_gh = DendDict['dend.ghbar_gh'] 
    g_pas = DendDict['dend.g_pas']
    e_pas = DendDict['dend.e_pas'] 
    htau_gh = GenDict['htau_gh']
    half_gh = DendDict['dend.half_gh']
    theta_m_L_Ca_inact = GenDict['theta_m_L_Ca_inact']
    tau_m_L_Ca_inact = GenDict['tau_m_L_Ca_inact']
    theta_h_L_Ca_inact = GenDict['theta_h_L_Ca_inact']
    tau_h_L_Ca_inact = GenDict['tau_h_L_Ca_inact']
    kappa_h_L_Ca_inact = GenDict['kappa_h_L_Ca_inact']

    return diam, L, Ra, cm, gcabar_L_Ca_inact, ghbar_gh, g_pas, e_pas, htau_gh, half_gh, theta_m_L_Ca_inact, \
tau_m_L_Ca_inact, theta_h_L_Ca_inact, tau_h_L_Ca_inact, kappa_h_L_Ca_inact

def Set_Parameters_D2(self, no_mu, filepath, CaPIC, i):
    
    SomaDict, DendDict, GenDict = Load_Parameters(self, no_mu, filepath, i)
        
    ## Set Biophysical Parameters   
    diam = DendDict['dend.diam']  
    L = DendDict['dend.L']
    Ra = DendDict['dend.Ra']
    cm = DendDict['dend.cm']
    gcabar_L_Ca_inact = CaPIC*(DendDict['d2.gcabar_L_Ca_inact'])
    ghbar_gh = DendDict['dend.ghbar_gh'] 
    g_pas = DendDict['dend.g_pas']
    e_pas = DendDict['dend.e_pas'] 
    htau_gh = GenDict['htau_gh']
    half_gh = DendDict['dend.half_gh']
    theta_m_L_Ca_inact = GenDict['theta_m_L_Ca_inact']
    tau_m_L_Ca_inact = GenDict['tau_m_L_Ca_inact']
    theta_h_L_Ca_inact = GenDict['theta_h_L_Ca_inact']
    tau_h_L_Ca_inact = GenDict['tau_h_L_Ca_inact']
    kappa_h_L_Ca_inact = GenDict['kappa_h_L_Ca_inact']

    return diam, L, Ra, cm, gcabar_L_Ca_inact, ghbar_gh, g_pas, e_pas, htau_gh, half_gh, theta_m_L_Ca_inact, \
tau_m_L_Ca_inact, theta_h_L_Ca_inact, tau_h_L_Ca_inact, kappa_h_L_Ca_inact

def Set_Parameters_D3(self, no_mu, filepath, CaPIC, i):
    
    SomaDict, DendDict, GenDict = Load_Parameters(self, no_mu, filepath, i)
        
    ## Set Biophysical Parameters   
    diam = DendDict['dend.diam']  
    L = DendDict['dend.L']
    Ra = DendDict['dend.Ra']
    cm = DendDict['dend.cm']
    gcabar_L_Ca_inact = CaPIC*(DendDict['d3.gcabar_L_Ca_inact'])
    ghbar_gh = DendDict['dend.ghbar_gh'] 
    g_pas = DendDict['dend.g_pas']
    e_pas = DendDict['dend.e_pas'] 
    htau_gh = GenDict['htau_gh']
    half_gh = DendDict['dend.half_gh']
    theta_m_L_Ca_inact = GenDict['theta_m_L_Ca_inact']
    tau_m_L_Ca_inact = GenDict['tau_m_L_Ca_inact']
    theta_h_L_Ca_inact = GenDict['theta_h_L_Ca_inact']
    tau_h_L_Ca_inact = GenDict['tau_h_L_Ca_inact']
    kappa_h_L_Ca_inact = GenDict['kappa_h_L_Ca_inact']

    return diam, L, Ra, cm, gcabar_L_Ca_inact, ghbar_gh, g_pas, e_pas, htau_gh, half_gh, theta_m_L_Ca_inact, \
tau_m_L_Ca_inact, theta_h_L_Ca_inact, tau_h_L_Ca_inact, kappa_h_L_Ca_inact

def Set_Parameters_D4(self, no_mu, filepath, CaPIC, i):
    
    SomaDict, DendDict, GenDict = Load_Parameters(self, no_mu, filepath, i)
        
    ## Set Biophysical Parameters   
    diam = DendDict['dend.diam']  
    L = DendDict['dend.L']
    Ra = DendDict['dend.Ra']
    cm = DendDict['dend.cm']
    gcabar_L_Ca_inact = CaPIC*(DendDict['d4.gcabar_L_Ca_inact'])
    ghbar_gh = DendDict['dend.ghbar_gh'] 
    g_pas = DendDict['dend.g_pas']
    e_pas = DendDict['dend.e_pas'] 
    htau_gh = GenDict['htau_gh']
    half_gh = DendDict['dend.half_gh']
    theta_m_L_Ca_inact = GenDict['theta_m_L_Ca_inact']
    tau_m_L_Ca_inact = GenDict['tau_m_L_Ca_inact']
    theta_h_L_Ca_inact = GenDict['theta_h_L_Ca_inact']
    tau_h_L_Ca_inact = GenDict['tau_h_L_Ca_inact']
    kappa_h_L_Ca_inact = GenDict['kappa_h_L_Ca_inact']

    return diam, L, Ra, cm, gcabar_L_Ca_inact, ghbar_gh, g_pas, e_pas, htau_gh, half_gh, theta_m_L_Ca_inact, \
tau_m_L_Ca_inact, theta_h_L_Ca_inact, tau_h_L_Ca_inact, kappa_h_L_Ca_inact        