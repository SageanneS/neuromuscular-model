# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 12:37:25 2019

@author: Sageanne Senneff
"""

import numpy as np
from random import uniform
           
def Excitatory_Cell_Generation(self, k, fr2, dt, no_mu, total_no_input_ind, durn, max_syn, a_s, tmax_II, tau_II, ts, d_n, EPSP_dt, start, endtime):
    
    ## EPSP Generation for Independent Excitatory Inputs
    t_x2 = np.arange(0, tmax_II + 1, EPSP_dt)
    gal2 = np.zeros([np.size(t_x2)])
    galp2 = tau_II/np.exp(1)
    tr2 = t_x2[round(ts/EPSP_dt)-1:len(t_x2)]-(ts-EPSP_dt)
    gal2[round(ts/EPSP_dt)-1:len(t_x2)] = (tr2 + d_n)*np.exp((-(tr2+d_n))/tau_II)/galp2
    EPSP = gal2
    
    pulse_train_ind = np.zeros([no_mu, total_no_input_ind, len(np.arange(0, durn, dt))])
    pulse_train_ind_ip = np.zeros([no_mu, total_no_input_ind, len(np.arange(0, durn, dt)) + len(EPSP)])
    pulse_count_ind = np.zeros([no_mu, total_no_input_ind]) 
    temp_pulse_ex = np.zeros([no_mu, len(np.arange(0, durn, dt))])
    
    ## Generate a 30 second signal
    temp1 = 0.1 + np.random.uniform(low = 0.1, high = 0.75, size=(total_no_input_ind, max_syn))
    for t in range(int(start), int(endtime)):
        for i in range(no_mu):
            for no_input in range(total_no_input_ind):    
                pspike = fr2*dt*k
                if (pspike > uniform(0,1)):
                    pulse_train_ind[i,no_input, t] = 1           # Fire with poisson distribution, renewal process
                else:
                    pulse_train_ind[i,no_input, t] = 0           # Fire with poisson distribution, renewal process
                if pulse_train_ind[i,no_input,t] == 1:            # Remove every kth Firing Time
                    pulse_count_ind[i,no_input] = pulse_count_ind[i,no_input] + 1
                    if np.remainder(pulse_count_ind[i,no_input], k) == 0:
                        pulse_train_ind[i,no_input,t] = 0
                    else:
                        pulse_train_ind[i,no_input,t] = sum(temp1[no_input,:] > np.random.uniform(low = 0, high = 1, size = max_syn), 1)*a_s
                        EPSPconj = EPSP.conj().transpose()
                        pulse_train_ind_ip[i,no_input, t:t + len(EPSP)] = np.squeeze(pulse_train_ind_ip[i,no_input, t:t + len(EPSP)]) + (EPSPconj*pulse_train_ind[i,no_input,t])            
            temp_pulse_ex[i,t] = sum(pulse_train_ind_ip[i,:,t]) 
            
    return temp_pulse_ex

def Inhibitory_Cell_Generation(self, k, fr2, dt, no_mu, total_no_input_inhib, durn, max_syn, a_s_inhib, tmax_II, tau_II, ts, d_n, EPSP_dt, start, endtime):
    
    ## EPSP Generation for Independent Inhibitory Inputs
    t_x2 = np.arange(0, tmax_II + 1, EPSP_dt)
    gal2 = np.zeros([np.size(t_x2)])
    galp2 = tau_II/np.exp(1)
    tr2 = t_x2[round(ts/EPSP_dt)-1:len(t_x2)]-(ts-EPSP_dt)
    gal2[round(ts/EPSP_dt)-1:len(t_x2)] = (tr2 + d_n)*np.exp((-(tr2+d_n))/tau_II)/galp2
    EPSP = gal2
    
    pulse_train_inhib = np.zeros([no_mu, total_no_input_inhib, len(np.arange(0, durn, dt))])
    pulse_train_inhib_ip = np.zeros([no_mu, total_no_input_inhib, len(np.arange(0, durn, dt)) + len(EPSP)])
    pulse_count_inhib = np.zeros([no_mu, total_no_input_inhib]) 
    temp_pulse_inhib = np.zeros([no_mu, len(np.arange(0, durn, dt))])
    
    ## Generate a 30 second signal
    temp2 = 0.1 + np.random.uniform(low = 0.1, high = 0.75, size=(total_no_input_inhib, max_syn))
    for t in range(int(start), int(endtime)):
        for i in range(self.no_mu):
            for no_input in range(total_no_input_inhib):    
                pspike = fr2*dt*k
                if (pspike > uniform(0,1)):
                    pulse_train_inhib[i,no_input, t] = 1                 # Fire with poisson distribution, renewal process
                else:
                    pulse_train_inhib[i,no_input, t] = 0                 # Fire with poisson distribution, renewal process
                if pulse_train_inhib[i,no_input,t] == 1:                  # Remove every kth Firing Time
                    pulse_count_inhib[i,no_input] = pulse_count_inhib[i,no_input] + 1
                    if np.remainder(pulse_count_inhib[i,no_input], k) == 0:
                        pulse_train_inhib[i,no_input,t] = 0
                    else:
                        pulse_train_inhib[i,no_input,t] = sum(temp2[no_input,:] > np.random.uniform(low = 0, high = 1, size = max_syn), 1)*a_s_inhib
                        EPSPconj = EPSP.conj().transpose()
                        pulse_train_inhib_ip[i,no_input, t:t + len(EPSP)] = np.squeeze(pulse_train_inhib_ip[i,no_input, t:t + len(EPSP)]) + (EPSPconj*pulse_train_inhib[i,no_input, t])            
            temp_pulse_inhib[i,t] = sum(pulse_train_inhib_ip[i,:,t])    

    return temp_pulse_inhib
