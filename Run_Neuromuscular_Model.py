# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 20:05:32 2019

@author: Sageanne Senneff
"""
from neuron import h
from Cortical_Cell_Generation import Excitatory_Cell_Generation, Inhibitory_Cell_Generation
from Gather_Parameters import Set_Parameters_Soma, Set_Parameters_D1, Set_Parameters_D2, Set_Parameters_D3, Set_Parameters_D4
from Correlated_Cortical_Cell_Generation import Correlated_Cell_Generation
from Force_Controller import Force_Feedback
import numpy as np
from scipy.stats import norm
import scipy.io as sio

class Run_Neuromuscular_Model():
    
    ## Load Necessary Neuron Libraries                             
    h.load_file("stdrun.hoc")                                                     
    h.load_file("stdlib.hoc")    
    #h('nrn_load_dll("C:/Users/Sageanne Senneff/Desktop/Code/Powers Model/nrnmech.dll")')

    def __init__(self):
        
        
        h.celsius = 37.0
        h.tstop = 5000.0
        h.dt = 0.025
        
        self.no_mu = 5
        self.dt = 0.001
        self.durn = 26.5
        self.CaPIC = 1.0
        self.filepath = 'C:/Users/Sageanne Senneff/Desktop/Code/Powers Model'
        
        self.pk = 0.2
        self.fr = 2.0
        self.fr2 = 6.2
        self.delaystart = 0.5
        self.ss_length = 20.0
        self.k = 2.0
        self.total_no_input_ind = 75
        self.total_no_input_inhib = 25
        self.max_syn = 10
        self.a_s = 68e-10                                                           
        self.a_s_inhib = (self.a_s/6)*5                                                      
        self.tmax_CI = 25
        self.tmax_II = 20
        self.tau_CI = 3
        self.tau_II = 2
        self.ts = 2
        self.d_n = 0.58
        self.EPSP_dt = 1
        self.start = 0
        self.endtime = self.durn/self.dt

        self.V_e = -65.0
        self.V_i = 0
        self.Fmax = 100000
        self.cortical_sig_amp = 0.02

    def __create__(self):

        self.all_soma = list()
        for i in range(self.no_mu):
            soma = h.Section()
            soma.insert('na3rp')
            soma.insert('naps')
            soma.insert('kdrRL')
            soma.insert('mAHP')
            soma.insert('gh') 
            soma.insert('pas') 
            
#            soma.diam, soma.L, soma.Ra, soma.cm, soma.gbar_na3rp, soma.gbar_naps, soma.gMax_kdrRL, soma.gcamax_mAHP, soma.gkcamax_mAHP, soma.ghbar_gh, \
#            soma.g_pas, soma.ek, soma.e_pas, soma.taur_mAHP, soma.mtauca_mAHP, soma.htau_gh, soma.tmin_kdrRL, soma.taumax_kdrRL, soma.sh_na3rp, \
#            soma.sh_naps,soma.ar_na3rp, soma.ar_naps, soma.half_gh, soma.qinf_na3rp, soma.thinf_na3rp, soma.vslope_naps, soma.asvh_naps, soma.bsvh_naps, \
#            soma.mvhalfca_mAHP, soma.mVh_kdrRL = Set_Parameters_Soma(self, self.no_mu, self.filepath, self.CaPIC, i) 
                        
            self.all_soma.append(soma)          
        
        self.all_d1 = list()
        for i in range(self.no_mu): 
            
            d1 = h.Section()
            d1.insert('pas')  
            d1.insert('L_Ca_inact')
            d1.insert('gh')
            d1.connect(soma(1),0)
            
            d1.diam, d1.L, d1.Ra, d1.cm, d1.gcabar_L_Ca_inact, d1.ghbar_gh, d1.g_pas, d1.e_pas, d1.htau_gh, d1.half_gh, d1.theta_m_L_Ca_inact, \
            d1.tau_m_L_Ca_inact, d1.theta_h_L_Ca_inact, d1.tau_h_L_Ca_inact, d1.kappa_h_L_Ca_inact = Set_Parameters_D1(self, self.no_mu, self.filepath, self.CaPIC, i) 
            
            self.all_d1.append(d1)    

        self.all_d2 = list()
        for i in range(self.no_mu): 
            d2 = h.Section()
            d2.insert('pas')  
            d2.insert('L_Ca_inact')
            d2.insert('gh')
            d2.connect(soma(1),0)
            
            d2.diam, d2.L, d2.Ra, d2.cm, d2.gcabar_L_Ca_inact, d2.ghbar_gh, d2.g_pas, d2.e_pas, d2.htau_gh, d2.half_gh, d2.theta_m_L_Ca_inact, \
            d2.tau_m_L_Ca_inact, d2.theta_h_L_Ca_inact, d2.tau_h_L_Ca_inact, d2.kappa_h_L_Ca_inact = Set_Parameters_D2(self, self.no_mu, self.filepath, self.CaPIC, i) 
            
            self.all_d2.append(d2)    

        self.all_d3 = list()    
        for i in range(self.no_mu): 
            d3 = h.Section()
            d3.insert('pas')  
            d3.insert('L_Ca_inact')
            d3.insert('gh')
            d3.connect(soma(0),0)
            
            d3.diam, d3.L, d3.Ra, d3.cm, d3.gcabar_L_Ca_inact, d3.ghbar_gh, d3.g_pas, d3.e_pas, d3.htau_gh, d3.half_gh, d3.theta_m_L_Ca_inact, \
            d3.tau_m_L_Ca_inact, d3.theta_h_L_Ca_inact, d3.tau_h_L_Ca_inact, d3.kappa_h_L_Ca_inact = Set_Parameters_D3(self, self.no_mu, self.filepath, self.CaPIC, i) 
            
            self.all_d3.append(d3)    

        self.all_d4 = list()
        for i in range(self.no_mu): 
            d4 = h.Section()
            d4.insert('pas')  
            d4.insert('L_Ca_inact')
            d4.insert('gh')
            d4.connect(soma(0),0)    
            
            d4.diam, d4.L, d4.Ra, d4.cm, d4.gcabar_L_Ca_inact, d4.ghbar_gh, d4.g_pas, d4.e_pas, d4.htau_gh, d4.half_gh, d4.theta_m_L_Ca_inact, \
            d4.tau_m_L_Ca_inact, d4.theta_h_L_Ca_inact, d4.tau_h_L_Ca_inact, d4.kappa_h_L_Ca_inact = Set_Parameters_D4(self, self.no_mu, self.filepath, self.CaPIC, i) 
            
            self.all_d4.append(d4) 

    def integrate(self):
        
        spike_times = []
        ipi = []
        self.Force = []
        
        while h.t < h.tstop:
            
            #for i in range(self.no_mu):
            if self.all_vsoma[0].x[int(h.t/h.dt)] >= - 20.0 and self.all_vsoma[0].x[int(h.t/h.dt)-1] < - 20.0:
                spike_times.append(h.t)
                if len(spike_times) > 2:
                    ipi.append(spike_times[-1]-spike_times[-2])
            
            for i in range(len(self.force_controller_call_times)-1):
                if np.round(h.t, decimals = 2) == self.force_controller_call_times[i]:
                    for j in range(len(ipi)):    
                        sig = 1 - np.exp(-2*((self.T/ipi[j])**3))       
                        if all(jj > 0.4 for jj in (self.T/ipi[j])):
                            g = sig/(self.T/ipi[j])
                            g = g/self.norm_twitch
                        else:
                            g = 1  
                        f = ((((g*self.P/self.T)*np.arange(0, 0.5, 0.005)))*np.exp(1 - (np.arange(0, 0.5, 0.005)/self.T)))
                        self.Force.append(f)  
                    
 
#            for j in range(len(self.force_controller_call_times)-1):
#                if np.round(h.t, decimals = 2) == self.force_controller_call_times[j]:
#                 sig = 1 - np.exp(-2*((T/ipi[j])**3))       
#                            if all(j > 0.4 for j in T/ipi[j]):
#                                g = sig/(T/ipi[j])
#                                g = g/norm_twitchs
#                            else:
#                                g = 1  
#                            f = ((((g*P/T)*np.arange(0, 0.5, 0.005)))*np.exp(1 - (np.arange(0, 0.5, 0.005)/T)))
#                            Force.append(f)           
 


                    #spike_times.append(h.t)
#                    spike_times = h.t
#                    self.all_spike_times.append(spike_times)
#                    if len(self.all_spike_times) > 2:
#                        ipi = self.all_spike_times[-1]-self.all_spike_times[-2]
#                        print(ipi)
                    
                   # self.all_spike_times.append(spike_times)
#                    if h.vsoma[i][t-1] <= -20 and h.vsoma[i][t] > -20: 
#                         st.append([h.trec[i][t]])
#                         if len(st) > 2:
                                
#                for j in range(len(self.all_spike_times[i])):
#                    ipi = self.all_spike_times[i] - self.all_spike_times[i][-1]
#                    self.all_ipi.append(ipi)
                                    
#            for j in range(len(self.force_controller_call_times)-1):
#                if np.round(h.t, decimals = 2) == self.force_controller_call_times[j]:
#                    print(h.t)
#                    print(self.all_ipi[:])

                    
            h.fadvance()
          
           # print(self.all_spike_times[0])
#            for j in range(len(self.force_controller_call_times)-1):
#                if np.round(h.t, decimals = 2) == self.force_controller_call_times[j]:
#                    print(h.t)

#                    Force = []
#                    for ii in range(self.no_mu):
#                         
#                        for jj in range(len(save_voltage[i,:])):
#                            ipi[jj] = self.all_t[ii][jj+1] - self.all_t[i][j]     
#                            sig = 1 - np.exp(-2*((T/ipi[j])**3))       
#                            if all(j > 0.4 for j in T/ipi[j]):
#                                g = sig/(T/ipi[j])
#                                g = g/norm_twitchs
#                            else:
#                                g = 1  
#                            f = ((((g*P/T)*np.arange(0, 0.5, 0.005)))*np.exp(1 - (np.arange(0, 0.5, 0.005)/T)))
#                            Force.append(f)           
#                   
#                    if self.force_controller_call_times[j] <= 500:
#                        for t in range(int(self.force_controller_call_times[j])):
#                            self.All_Force[i, t] = np.mean(Force)
#                    else:
#                        for t in range(int(self.force_controller_call_times[j-1]), int(self.force_controller_call_times[j])):
#                            self.All_Force[i, int(self.force_controller_call_times[j-1]):int(self.force_controller_call_times[j])] = np.mean(Force)
#                    
#                    
#                    if self.pk == 1.0:
#                        self.fr = max_fr
#                        self.fr2 = max_fr2
#                    else:
#                        if self.force_controller_call_times[j] <= (cornerb-0.5)/self.dt:
#                           self.fr = 2.0
#                           self.fr2 = 6.2
#                        else:
#                            if self.force_controller_call_times[j]  >= (cornerb + adj_force + 2)/self.dt and self.force_controller_call_times[j]  < (cornerd-adj_force)/self.dt:
#                                err = (np.mean(self.All_Force[:, int(self.force_controller_call_times[j]) -(adj_force/self.dt):int(self.force_controller_call_times[j])]) - np.mean(target[:,int(int(self.force_controller_call_times[j])  + self.delaystart/self.dt + (adj_int/self.dt)/2)]))*(self.Fmax/np.sum(self.All_Force[:, int(self.force_controller_call_times[j])])) 
#                            elif self.force_controller_call_times[j]  > (cornerd-adj_force)/self.dt:
#                                err = np.sum(self.All_Force[:, int(self.force_controller_call_times[j])])/np.sum(self.All_Force[:, int(self.force_controller_call_times[j])])
#                            else:
#                                err=(np.sum(self.All_Force[:, int(self.force_controller_call_times[j])])-(self.pk*self.Fmax))/np.sum(self.All_Force[:, int(self.force_controller_call_times[j])])
#                            
#                            errmod1 = 0.015*(1/np.sqrt(self.pk))
#                            errmod2 = 0.2             
#                            fr2increase = 1.5
#                                                    
#                            # Check the force and increase firing rates if necessary
#                            if np.sum(self.All_Force[:, int(self.force_controller_call_times[j])])/self.Fmax <= self.pk-(errmod1*self.pk) and self.force_controller_call_times[j] < (cornerb+2)/self.dt:
#                                if (self.fr-self.fr*(errmod1*err)) > max_fr:
#                                    self.fr  = max_fr
#                                    self.fr2 = max_fr2
#                                else:
#                                    self.fr  = self.fr-self.fr*(errmod1*err)
#                                    self.fr2 = self.fr2-self.fr2*(errmod1*err*fr2increase)
#                            elif np.equal(self.force_controller_call_times[j]%(comp_int/self.dt), 0):
#                                if (self.fr-self.fr*(errmod2*err)) > max_fr:
#                                    self.fr  = max_fr
#                                    self.fr2 = max_fr2
#                                else:
#                                    self.fr  = self.fr-self.fr*(errmod2*err)
#                                    self.fr2 = self.fr2-self.fr2*(errmod2*err*fr2increase)         
#                            if  self.fr2 > max_fr2:
#                                self.fr2 = max_fr2         
#                        
#                    # Calculate new Cortical Signals         
#                    self.new_excitatory_cells = Excitatory_Cell_Generation(self, k = 2.0, fr2 = self.fr2, dt = self.dt, no_mu = self.no_mu, total_no_input_ind = 75, \
#                                                                     durn = self.durn, max_syn = 10, a_s = 68e-10, tmax_II = 20, tau_II = 2, ts = 2, d_n = 0.58, EPSP_dt = 1, \
#                                                                     start = self.force_controller_call_times[j], endtime = self.force_controller_call_times[j+1])
#                    self.new_inhibitory_cells = Inhibitory_Cell_Generation(self, k = 2.0, fr2 = self.fr2, dt = self.dt, no_mu = self.no_mu, total_no_input_inhib = 25, \
#                                                                     durn = self.durn, max_syn = 10, a_s_inhib = (68e-10/6)*5, tmax_II = 20, tau_II = 2, \
#                                                                     ts = 2, d_n = 0.58, EPSP_dt = 1, start = self.force_controller_call_times[j], \
#                                                                     endtime = self.force_controller_call_times[j+1])
#                    self.new_correlated_cells = Correlated_Cell_Generation(self, no_mu = self.no_mu, filepath = self.filepath, durn = self.durn, dt = self.dt, fr = self.fr, fr2 = self.fr2, \
#                                                                     tmax_CI = 25, tau_CI = 3, ts = 2, d_n = 0.58, EPSP_dt = 1, cortical_sig_amp = 0.01, pk = self.pk, tmax_II = 20, tau_II = 2, \
#                                                                     start = self.force_controller_call_times[j], endtime = self.force_controller_call_times[j+1])
#                    # Update MN Inputs
#                    for i in range(self.no_mu):
#                       for ii in [self.all_d1[i], self.all_d2[i], self.all_d3[i], self.all_d4[i]]:
#                            updated_cortical_signals[int(self.force_controller_call_times[j]):int(self.force_controller_call_times[j+1])] = \
#                            self.new_excitatory_cells[i][int(self.force_controller_call_times[j]):int(self.force_controller_call_times[j+1])]*(ii(0.5).v - self.V_e)\
#                            + self.new_inhibitory_cells[i][int(self.force_controller_call_times[j]):int(self.force_controller_call_times[j+1])]*(ii(0.5).v - self.V_i)\
#                            + self.new_correlated_cells[i][int(self.force_controller_call_times[j]):int(self.force_controller_call_times[j+1])]*(ii(0.5).v - self.V_e)

    def go(self):
        
        h.finitialize(-65)
        self.integrate()
        
        
    def __main__(self):    
        
        ## Create Sections and Set Parameters for MUs
        self.__create__()

        ## Set Up Cortical Signals
        self.all_excitatory = Excitatory_Cell_Generation(self, self.k, self.fr2, self.dt, self.no_mu, self.total_no_input_ind, \
                     self.durn, self.max_syn, self.a_s, self.tmax_II, self.tau_II, self.ts, self.d_n, self.EPSP_dt, self.start, self.endtime)    
        self.all_inhibitory = Inhibitory_Cell_Generation(self, self.k, self.fr2, self.dt, self.no_mu, self.total_no_input_inhib, \
                        self.durn, self.max_syn, self.a_s_inhib, self.tmax_II, self.tau_II, self.ts, self.d_n, self.EPSP_dt, self.start, self.endtime) 
        self.all_correlated = Correlated_Cell_Generation(self, self.no_mu, self.filepath, self.durn, self.dt, self.fr, self.fr2, self.tmax_CI, self.tau_CI, \
                     self.ts, self.d_n, self.EPSP_dt, self.cortical_sig_amp, self.pk, self.tmax_II, self.tau_II, self.start, self.endtime)
        
        h('objref iclamp, cortical_signals, tvec')
        for i in range(self.no_mu):
            for ii in [self.all_d1[i], self.all_d2[i], self.all_d3[i], self.all_d4[i]]:
                h.cortical_signals = h.Vector(h.tstop)
                h.tvec = h.Vector(h.tstop)
                for t in range(int(h.tstop)):
                    h.cortical_signals.x[t] = 100 #self.all_excitatory[i,t] + self.all_inhibitory[i,t] + self.all_correlated[i,t]
                    h.tvec.x[t] = t
                h.iclamp = h.IClamp(ii(0.5))
                h.iclamp.dur = h.tstop
                h.iclamp.delay = 0
                h.cortical_signals.play(h.iclamp, h.iclamp._ref_amp, h.tvec, 1)

        self.All_Force = np.zeros([self.no_mu, len(np.arange(0, self.durn, self.dt))]) 
        
        self.all_time = list()
        for i in range(self.no_mu):
            time = h.Vector()
            time.record(h._ref_t)
            self.all_time.append(time)
    
        self.all_vsoma = list()
        for i in range(self.no_mu):
            vsoma = h.Vector()
            vsoma.record(self.all_soma[i](0.5)._ref_v)
            self.all_vsoma.append(vsoma)
            
        adj_int = 0.15
        adj_force = 1.0 
        max_fr = 3.2                                                           # Max Cortical Firing Rate
        max_fr2 = 50     
        comp_int = 0.75 

        ## Set up Force Trajectory
        cornera = 0 + self.delaystart                                         
        cornerb = (self.pk*10) + self.delaystart
        cornerc = (self.pk*10) + self.ss_length + self.delaystart
        cornerd = (self.pk*10*2) + self.ss_length +self. delaystart
        target = np.zeros([self.no_mu, len(np.arange(0, self.durn, self.dt))])
        target[:,int(cornera/self.dt):int(cornerb/self.dt)] = np.arange(self.pk/((cornerb-cornera)/self.dt), self.pk, self.pk/((cornerb-cornera)/self.dt))
        target[:,int(cornerb/self.dt):int(cornerc/self.dt)] = np.ones([1, len(np.arange(cornerb/self.dt, cornerc/self.dt, 1))])*self.pk   
        target[:,int(cornerc/self.dt):int(cornerd/self.dt)] = np.arange(self.pk, 0, -self.pk/((cornerd-cornerc)/self.dt)) 

        self.P = np.zeros([100, 1])
        self.T = np.zeros([100, 1])
        timedist = np.arange(0,1 + 0.001, 1/(100 + 1))
        pd = np.random.weibull(1, 102)
        Tdistrib = norm.cdf(pd,timedist)
        Tdistrib = np.array(Tdistrib).tolist()
        Tdistrib.pop()
        Tdistrib.pop(0)
        Tdistrib = np.array(Tdistrib)
        Tdistrib = Tdistrib*40 + 31
        timedist = np.arange(0,5 + 0.001, 5/100)
        Pdistrib = (np.exp(timedist)/150)*70 + 1
        for m in range(100):
            self.T[m] = Tdistrib[m]   
            self.P[m] = Pdistrib[m]
        self.T = np.fliplr(self.T*1e-3)       
        self.norm_twitch = 1 - np.exp(-2*((0.4)**3))

        self.force_controller_call_times = np.arange(self.delaystart*1000, h.tstop, 150.0)                
        self.all_spike_times = list()
        self.all_ipi = list()

        self.go()
        
        print(self.all_spike_times)
        
        ## Plot the Data
        from matplotlib import pyplot as plt
        plt.figure(0)
        plt.plot(np.arange(0, h.tstop+h.dt+h.dt, h.dt), self.all_vsoma[0], 'k') 
        plt.figure(1)
        plt.plot(self.Force, 'k') 
                  
        myData = sio.savemat('Model_Output.mat', {'vsoma':self.all_vsoma, 'spike_times':self.all_spike_times, 'IPI': self.all_ipi})
        
        return myData
        
if __name__ == '__main__':
    runner = Run_Neuromuscular_Model()
    runner.__main__()
                
