# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 19:00:14 2020

@author: kevin
"""

#this code is to examine the evolution of weights for different information sources: outcome vs. action.

import matplotlib.pyplot as plt
import random
import numpy as np

N=5000
T1_percent=0.9
s1=1 #s represent the technological efficacy of T1 and T2. Larger s means higher fitness.
s2=2
beta=0 #represents the strength of genetic selection. beta=0 means all individuals have the same fitness
generation=300
mut_mag=0.1
sample_size=10
ini_w_o=1
D=0.6

gamma=0.5

fit_function="normal"
prob_method="bayesian" #can be linear or bayesian


def fitness(agent):
    if fit_function=="normal":
        if agent.T=="T1":
            return s1**beta
        elif agent.T=="T2":
            return s2**beta
    if fit_function=="extreme":
        if agent.T=="T1":
            return 0
        if agent.T=="T2":
            return 1
def mutate(agent):
    #agent.w_o=agent.w_o+np.random.normal(0,mut_mag)
    agent.w_a=agent.w_a+np.random.normal(0,mut_mag)
    
def T1_freq(pop):
    T1_num=0
    for i in pop:
        if i.T=="T1":
            T1_num=T1_num+1
    return T1_num/len(pop)

def ave_weight(pop,w):
    w_o_list=[]
    w_a_list=[]
    if w=="w_o":
        for i in pop:
            w_o_list.append(i.w_o)
        return np.mean(w_o_list)
    if w=="w_a":
        for i in pop:
            w_a_list.append(i.w_a)
        return np.mean(w_a_list)

class agent:
    def __init__(self,T,w_o,w_a):
        self.T=T
        self.w_o=w_o
        self.w_a=w_a


#set up initial population

def main(T1_percent,s1,s2,ini_w_o,sample_size,D):
    pop=[]
    
    for i in range(N):
        if i<N*T1_percent:#fix w_o, only let w_a evolve
            i=agent(T="T1",w_o=ini_w_o,w_a=np.random.normal(1,0.1))
        else:
            i=agent(T="T2",w_o=ini_w_o,w_a=np.random.normal(1,0.1))
        pop.append(i)
    
    
    timer=0
    global T1_freq_time_series; T1_freq_time_series=[]
    global ave_w_a_time_series; ave_w_a_time_series=[]
    
    
    for gen in range(generation):
        T1_freq_time_series.append(T1_freq(pop))
        ave_w_a_time_series.append(ave_weight(pop, "w_a"))
        fit_list=[] #create a fitness list of all the individuals in the 
        for ind in pop:
            fit_list.append(fitness(ind))
        fit_list_normalized=[]
        for i in fit_list:
            fit_list_normalized.append(i/sum(fit_list))
        #then create a offspring_pop:
        pop_f1=[]
        #pick individuals in pop to reproduce
        for i in range(N):
            parent_indx=np.random.choice(range(N),size=None,p=fit_list_normalized)
            pop_f1.append(agent(T=[],w_o=pop[parent_indx].w_o,w_a=pop[parent_indx].w_a)) #offspring inherit weights from parents
    
        #weights then mutate:
        for i in pop_f1:
            i=mutate(i)
    
        #now, agents in pop_f1 sample some number of individuals from pop to decide what trait to adopt:
        for i in pop_f1:
            model_list=[]
            for j in range(int(sample_size)):
                model_list.append(pop[random.randint(0,N-1)])
    
            T1_num=0
            T2_num=0
            for model in model_list:
                if model.T=="T1":
                    T1_num=T1_num+1
                elif model.T=="T2":
                    T2_num=T2_num+1
    
            if T1_num==0: #no T1 ind in model set
                p_T1_adopted=0
            elif T2_num==0:
                p_T1_adopted=1
                
            elif prob_method=="linear":
                if T1_num>sample_size/2:
                    p_T1_adopted=gamma*(T1_num+D)/sample_size+(1-gamma)*s1/(s1+s2)
                elif T1_num==sample_size/2:
                    p_T1_adopted=gamma*(T1_num)/sample_size+(1-gamma)*s1/(s1+s2)
                elif T1_num<sample_size/2:
                    p_T1_adopted=gamma*(T1_num-D)/sample_size+(1-gamma)*s1/(s1+s2)
            elif prob_method=="bayesian": #a mixture of T1 and T2 ind in model set
                if T1_num>sample_size/2:
                    p_T1_adopted=((T1_num+D)*i.w_a+s1*i.w_o)/(sample_size*i.w_a+(s1+s2)*i.w_o)
                elif T1_num==sample_size/2:
                    p_T1_adopted=((T1_num)*i.w_a+s1*i.w_o)/(sample_size*i.w_a+(s1+s2)*i.w_o)
                elif T1_num<sample_size/2:
                    p_T1_adopted=((T1_num-D)*i.w_a+s1*i.w_o)/(sample_size*i.w_a+(s1+s2)*i.w_o)
            
            if p_T1_adopted>random.random():
                i.T="T1"
            else:
                i.T="T2"
        pop=pop_f1
        timer=timer+1
        print ("generation=", timer)
    
    fig,ax1=plt.subplots()
    ax1.set_xlabel("generation")
    ax1.set_ylabel("T1 frequency")    
    plot1=ax1.plot(range(generation),T1_freq_time_series,color="b",label="T1 frequency")
    plt.axhline(y=s1/(s1+s2), color='black', linestyle='--')
    plt.ylim(0,1)
    
    ax2=ax1.twinx()
    ax2.set_ylabel("action weight")
    plot2=ax2.plot(range(generation),ave_w_a_time_series,color="r",label="action weight")
    plt.ylim(0,2)
    
    plots=plot1+plot2
    label_plots=[p.get_label() for p in plots]
    ax1.legend(plots,label_plots,loc=0)
    plt.savefig("D="+ str(D) +"N="+ str(N) + " T1_perc"+ str(T1_percent) + "s1="+ str(s1) + "s2="+ str(s2)+" samsize="+ str(sample_size) + ".png", format='png', dpi=1200)
    
    print ("D= ", D, "T1_percent= ", T1_percent)
                

D_set=[0,0.3,0.6]
D_array=[] #for T1_freq_time_series
for i in D_set:
    main(T1_percent, s1, s2, ini_w_o, sample_size,i)
    D_array.append(T1_freq_time_series)

plt.figure()
plt.style.use('ggplot')
plt.plot(range(generation),D_array[0],color="cyan",label="D=0.0")
plt.plot(range(generation),D_array[1],color="blue",label="D=0.3")
plt.plot(range(generation),D_array[2],color="darkblue",label="D=0.6")
plt.ylim(0,1)
plt.xlabel("Generation",fontsize=14)
plt.ylabel("Frequency",fontsize=14)
plt.axhline(y=s1/(s1+s2), color='black', linestyle='--')

#plt.legend()



