# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 19:47:48 2019

@author: Sageanne Senneff
"""

import numpy as np       

def Force_Feedback(self, Fmax, delaystart, ss_length, pk, no_mu, durn, dt, all_force_twitches, t):

    Force.append(sum(all_force_twitches))
    print(Force)
    
    adj_int = 0.15
    adj_force = 1.0 
    max_fr = 3.2                                                               # Max Cortical Firing Rate
    max_fr2 = 50     
    comp_int = 0.75 

    ## Setup Force Trajectory
    cornera = 0 + delaystart                                                   # time (seconds) when max reached
    cornerb = (pk*10) + delaystart
    cornerc = (pk*10) + ss_length + delaystart
    cornerd = (pk*10*2) + ss_length + delaystart
    target = np.zeros([no_mu, len(np.arange(0, durn, dt))])
    target[:,int(cornera/dt):int(cornerb/dt)] = np.arange(pk/((cornerb-cornera)/dt), pk, pk/((cornerb-cornera)/dt))
    target[:,int(cornerb/dt):int(cornerc/dt)] = np.ones([1, len(np.arange(cornerb/dt, cornerc/dt, 1))])*pk   
    target[:,int(cornerc/dt):int(cornerd/dt)] = np.arange(pk, 0, -pk/((cornerd-cornerc)/dt)) 
    
    t = int(t)
    print(t)
    
    if pk == 1.0:
        fr = max_fr
        fr2 = max_fr2
    else:
        if t <= (cornerb-0.5)/dt:
            fr = 2.0
            fr2 = 6.2
        else:
            if t >= (cornerb + adj_force + 2)/dt and t < (cornerd-adj_force)/dt:
                err = (np.mean(all_force_twitches) - np.mean(target[:,int(t + delaystart/dt + (adj_int/dt)/2)]))*(Fmax/Force) 
            elif t > (cornerd-adj_force)/dt:
                err = Force/Force 
            else:
                err=(Force-(pk*Fmax))/Force   
            
            errmod1 = 0.015*(1/np.sqrt(pk))
            errmod2 = 0.2             
            fr2increase = 1.5
                                    
            # Check the force and increase firing rates if necessary
            if Force/Fmax <= pk-(errmod1*pk) and t < (cornerb+2)/dt:
                if (fr-fr*(errmod1*err)) > max_fr:
                    fr  = max_fr
                    fr2 = max_fr2
                else:
                    fr  = fr-fr*(errmod1*err)
                    fr2 = fr2-fr2*(errmod1*err*fr2increase)
            elif np.equal(t%(comp_int/dt), 0):
                if (fr-fr*(errmod2*err)) > max_fr:
                    fr  = max_fr
                    fr2 = max_fr2
                else:
                    fr  = fr-fr*(errmod2*err)
                    fr2 = fr2-fr2*(errmod2*err*fr2increase)         
            if  fr2 > max_fr2:
                fr2 = max_fr2 
    
    return fr, fr2

