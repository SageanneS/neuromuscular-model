# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 18:11:08 2019

@author: Sageanne Senneff
"""
import numpy as np
import scipy.io as sio
import random

def Pulse_Encoder(self, fr, fr2, total_no_input, cortical_sig, durn, dt, tmax_CI, tau_CI, ts, d_n, EPSP_dt, tmax_II, tau_II, start, endtime):
    
    ## Pulse Encoder Parameters
    ZMAX    = 1.0                                                              # Threshold at which ENCODER fires
    GAIN    = 1.0                                                              # Gain of the ENCODER
    H       = dt                                                               # Sampling Rate
    TAU     = 0.025                                                            # ENCODER time constant (Reduced value to increase encoder synchrony)
    coeff01 = TAU/(H+TAU)
    coeff02 = GAIN*H/(H+TAU)

    ## EPSP Generation for Cortical Input
    t_x = np.arange(0, tmax_CI + 1, EPSP_dt)
    gal = np.zeros([np.size(t_x)])
    galp = tau_CI/np.exp(1)
    tr = t_x[round(ts/EPSP_dt)-1:len(t_x)] - (ts - EPSP_dt)
    gal[round(ts/EPSP_dt)-1:len(t_x)] = (tr + d_n)*np.exp((-(tr+d_n))/tau_CI)/galp
    EPSP_CT = gal
    
    ## EPSP Generation for Independent Excitatory Inputs
    t_x2 = np.arange(0, tmax_II + 1, EPSP_dt)
    gal2 = np.zeros([np.size(t_x2)])
    galp2 = tau_II/np.exp(1)
    tr2 = t_x2[round(ts/EPSP_dt)-1:len(t_x2)]-(ts-EPSP_dt)
    gal2[round(ts/EPSP_dt)-1:len(t_x2)] = (tr2 + d_n)*np.exp((-(tr2+d_n))/tau_II)/galp2
    EPSP = gal2
    
    znow = np.zeros([total_no_input,1])
    tnow = np.zeros([total_no_input,1]) 
    zold = np.zeros([total_no_input,1])
    start_pe = np.ones([total_no_input,1])
    nspike = np.ones([total_no_input,1])
    pulse_train = np.zeros([len(np.arange(0, durn, dt)), total_no_input])
    pulse_train_total = np.zeros([total_no_input, len(np.arange(0, durn, dt)) + len(EPSP_CT)])

    for t in range(int(start), int(endtime)):
        for no_input in range(total_no_input):
            if (znow[no_input,0] < ZMAX):                                          # ENCODER below threshold
                zold[no_input,0] = znow[no_input,0]
                start_pe = 1
                if start_pe:
                    v1 = 2.0*random.uniform(0,1)-1.0
                    v2 = 2.0*random.uniform(0,1)-1.0
                    w = v1*v1 + v2*v2 
                    w = np.log(w)/w
                    if np.sign(w) == -1:
                        w = np.sqrt(-w-w)
                    else:
                        w = np.sqrt(w+w)
                    g1 = v1*w
                    g2 = v2*w
                    start = 0
                    normal_output = fr + fr2*g1
                else:
                    start_pe = 1
                    normal_output = fr + fr2*g2                   
                znow[no_input,0] = coeff01*znow[no_input,0] + coeff02*normal_output + cortical_sig[0,t-1]
                tnow[no_input,0] = tnow[no_input,0] + H
                pulse_train[t,no_input] = 0
            elif (znow[no_input,0] == ZMAX):                                       # ENCODER on threshold
                tnow[no_input,0] = 0
                znow[no_input,0] = 0
                zold[no_input,0] = 0
                nspike[no_input,0] = nspike[no_input,0] + 1
                pulse_train[t,no_input] = 1
                pulse_train_total[no_input,t:t+len(EPSP)] = (pulse_train_total[no_input,t:t+len(EPSP)])+(EPSP*pulse_train[t,no_input])
            else:                                                                  # ENCODER above threshold
                tmp = H*(znow[no_input,0]-ZMAX)/(znow[no_input,0]-zold[no_input,0])
                nspike[no_input,0]=nspike[no_input,0]+1
                if znow[no_input,0] == 0:
                    znow[no_input,0] = 0
                else:
                    znow[no_input,0] = znow[no_input,0]%ZMAX + (ZMAX*(np.sign(znow[no_input,0]) - 1)/2)
                zold[no_input,0] = znow[no_input,0]
                tnow[no_input,0] = tmp
                pulse_train[t,no_input] = 1
                pulse_train_total[no_input,t:t+len(EPSP)] = (pulse_train_total[no_input,t:t+len(EPSP)])+(EPSP*pulse_train[t,no_input])  
    return pulse_train_total
        
def Correlated_Cell_Generation(self, no_mu, filepath, durn, dt, fr, fr2, tmax_CI, tau_CI, ts, d_n, EPSP_dt, cortical_sig_amp, pk, tmax_II, tau_II, start, endtime):

    ## Connectivity Parameters
    Connectivity_Dict = sio.loadmat(filepath + '/Cortical Model/InputConnectionDirect_15pc_600_100.mat')
    InputConnection = np.array(Connectivity_Dict['InputConnection'])
    total_no_input = np.amax(InputConnection)                                  # total number of correlated inputs to each MN 
    
    ## Initialize Cortical Signals
    Beta_Dict = sio.loadmat(filepath + '/Cortical Model/cortical_sig_2000_highBeta.mat')    
    cortical_sig_temp = np.array(Beta_Dict['cortical_sig_temp'])
    cortical_sig = np.zeros([total_no_input, np.size(cortical_sig_temp[1])])    
    NoiseRatioCT = 1
    if (pk==0.1):
        for ix in range(1,2):
            cortical_sig = cortical_sig_temp[ix,:]/np.std(cortical_sig_temp[ix,:])*cortical_sig_amp                                    
            cortical_sig = np.tile(cortical_sig, [total_no_input,1])
        for ix in range(0,np.size(cortical_sig,0)):
            cortical_sig[ix,:] = cortical_sig[ix,:] + np.random.uniform(low = 0, high = 1, size = np.size(cortical_sig_temp[1]))*np.std(cortical_sig[ix,:])*(NoiseRatioCT*0.75)
            cortical_sig[ix,:] = cortical_sig[ix,:]/np.std(cortical_sig[ix,:])*cortical_sig_amp
    else:                        
        for ix in range(1,2):
            cortical_sig = cortical_sig_temp[ix,:]/np.std(cortical_sig_temp[ix,:])*cortical_sig_amp
        cortical_sig = np.tile(cortical_sig, [total_no_input,1])  
        for ix in range(0,np.size(cortical_sig,0)): 
            cortical_sig[ix,:] = cortical_sig[ix,:] + np.random.uniform(low = 0, high = 1, size = np.size(cortical_sig_temp[1]))*np.std(cortical_sig[ix,:])*NoiseRatioCT
            cortical_sig[ix,:] = cortical_sig[ix,:]/np.std(cortical_sig[ix,:])*cortical_sig_amp  
 
    pulse_train_total = Pulse_Encoder(self, fr, fr2, total_no_input, cortical_sig, durn, dt, tmax_CI, tau_CI, ts, d_n, EPSP_dt, tmax_II, tau_II, start, endtime)
    
    temp_pulse_CT = np.zeros([no_mu, len(np.arange(0, durn, dt))])                 # Common Modulatory Pulse Train
    for t in range(int(start), int(endtime)):
        for i in range(self.no_mu):
            temp_pulse_CT[i,t] = sum(pulse_train_total[InputConnection[:, i]-1,t])
                       
    return temp_pulse_CT
