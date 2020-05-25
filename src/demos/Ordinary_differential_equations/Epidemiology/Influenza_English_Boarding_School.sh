#! /bin/sh
#  Author:  April 2020: John Stockwell
# copyright 2020  John Stockwell

# SIR model

# Important quantities:						
# r0 = number of new infections per single infected host  	
#  1 < r0 < 1.5 for influenza, (2.2 to 2.7 for Covid-19), 12 to	
# 18 for measles.						
#  b, k, S0, and r0 are related via				
#  k = bS0/r0 = b/r0 when S0/N and S0=N 			
#  								
#  It is often easier to determine the recovery rate k (in units
#  of h and to determine reasonable estimate of S0 and of r0 	
#  and to calculate the infection rate b = kr0/S0 or b=kr0	
#  when S0=N and is normalized by N.				
#								
# S = total number susceptible to the infection			
# I = total number of those capable of passing on the infection	
# R = total number removed = dead + recovered			
#								
# When xi is nonzero, then there is a potential that fraction of 
# the removed population can be reinfected.			

 
#  The SRI model describes an epidemic in terms of
#    S = susceptibles in a population
#    I = infectives in a population
#    R = removed = recovered + dead
# 
#    S0 = initial value of the susceptibles
#    I0 = initial value of the infectives
#    R0 = initial removed value = 0
#    
#    S(t) + I(t) + R(t) = S0 + I0   = N for the unnormalized case.
#    If normalized by total population N, then S(t) + I(t) + R(t) = 1 
#    and S(t) starts at its maxium value of S0/N.   
#    
#    r0 = bS0/k  = basic reproduction rate
#    b = rate of infection
#    k = rate removal = recovery rate + death rate
#    xi = re-infection rate 
#    mu = birth rate  
#    nu = death rate
#     
#    The encounters between susceptibles and the infectives is represented
#    by the product SI  
# 
#   SIR model:  
# 	S'(t) =  - bSI 
# 	I'(t) = bSI- kI 
# 	R'(t) = kI 
#     
#   SIR model with vital statistics (mu birth rate, nu death rate):  
# 	S'(t) = mu - nuS - bSI 
# 	I'(t) = bSI - kI - nuI 
# 	R'(t) = kI -  nuR
# 
#   SIRS model with vital statistics (mu birth rate, nu death rate) and reinfection:  
# 	S'(t) = mu - nuS + xiR - bSI 
# 	I'(t) = bSI - kI - nuI 
# 	R'(t) = kI - xiR - nuR
# 
#  S(t)= susceptible members 
#  I(t)= infectives
#  R(t)= removed members = recovered + dead + sequestered
# 
#  There is an impiled flow from S(t) -> I(t) -> R(t), though infected
#  who are quarantined immediately become part of R(t). 
# 
#  The product xiR are the reinfected members of the recovered group, and are thus 
#  removed from the recovered group and fed back to the susceptible group.
#  
#  The product bSI denotes the interaction of the infective population with
#  the susceptible population..
# 
								
### SIR model 								
# Influenza in an English boarding school, 1978:		
# N=762 I0=1,  2 students infected per day, 1/2 of the infected	
# population removed per day. Take b=2 k=0.5 			
#								
# Normalized by N:						
 sir_epidemic h=0.1 stepmax=200 I0=1 b=2 k=.5 N=762 mode=SIR |
 xgraph n=200 nplot=3 d1=.1 style=normal label1="days" \
 title="Influenza: English boarding school" \
 width=1000 height=1000 &

# real data

cat data.txt | awk '{ print $2 }' | a2b n1=1 > real_data.bin

				
# Normalized by N, output scaled by N:
 sir_epidemic h=.1 stepmax=200 I0=1 b=2 k=.5 N=762 mode=SIR scale=1 > model_data.bin

cat real_data.bin model_data.bin |
 xgraph n=15,200,200,200 nplot=4 d1=1,.1,.1,.1 linewidth=0,3,3,3 linecolor=4,5,4,6  \
 marksize=8,0,0,0 mark=3,0,0,0 \
 style=normal label1="days" \
 title="Model and Data: English boarding school" \
 width=1000 height=1000  \
 label2="number of students" &				

exit 0
