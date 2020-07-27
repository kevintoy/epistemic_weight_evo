# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 16:12:57 2020

@author: kevin
"""
import time
start_time = time.time()
import numpy as np
import random
import matplotlib.pyplot as plt

#set up a simulation to check analytic probability formulations
w_a=1
w_o=0
s1=1
s2=1
D=0.6
N=1000
sample_size=10
T1_percent=0.9

def prob_dist(T1_adopt_prob_list):
    plt.hist(T1_adopt_prob_list,bins=1000)

def T1_prop(pop): #find the proportion of types of cultural models in naive individuals' model set
    global pop_T1_freq_list; pop_T1_freq_list=[]
    for i in pop:
        pop_T1_freq_list.append(i.count("T1"))
    return pop_T1_freq_list

def T1_adopt_prob(pop_prop):
    global T1_adopt_prob_list; T1_adopt_prob_list=[]
    for i in pop_prop:
        if i==0:
            T1_adopt_prob_list.append(0)
        elif i==sample_size:
            T1_adopt_prob_list.append(1)
        elif i>sample_size/2:
            T1_adopt_prob_list.append(((i+D)*w_a+s1*w_o)/(sample_size*w_a+(s1+s2)*w_o))
        elif i<sample_size/2:
            T1_adopt_prob_list.append(((i-D)*w_a+s1*w_o)/(sample_size*w_a+(s1+s2)*w_o))
        elif i==sample_size/2:
            T1_adopt_prob_list.append(((i)*w_a+s1*w_o)/(sample_size*w_a+(s1+s2)*w_o))
    return T1_adopt_prob_list
            
def T1_adopt_decision(adopt_prob_list):
    
    global new_pop; new_pop=[]
    for i in adopt_prob_list:
        if i>random.random():
            new_pop.append("T1")
        else:
            new_pop.append("T2")
    return new_pop        

def new_T1_freq(new_pop):
    return new_pop.count("T1")/N
    print (new_pop.count("T1")/N)
    
def main(T1_percent):
    pop=[]
    
    for i in range(N):
        if i<N*T1_percent:#fix w_o, only let w_a evolve
            i="T1"
        else:
            i="T2"
        pop.append(i)
    
    pop_f1=[]
    
    for i in pop:
        model_list=[]
        for j in range(int(sample_size)):
            model_list.append(pop[random.randint(0,N-1)])
        pop_f1.append(model_list)
    
    
    #now check the T1 frequency in the next generation
    adoption_probability=T1_adopt_prob(T1_prop(pop_f1))
    
    return new_T1_freq(T1_adopt_decision(adoption_probability))
    
T1_freq_t1 = np.linspace(0,1,num=1001)
T1_freq_t2=[]
timer=0
for i in T1_freq_t1:
    T1_freq_t2.append(main(i))
    timer=timer+1
    print ("timer=",timer)

plt.figure()
plt.plot(T1_freq_t1,T1_freq_t2,linestyle="-")


print("--- %s seconds ---" % (time.time() - start_time))
