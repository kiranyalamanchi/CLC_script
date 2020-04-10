##Mechanism check for violation of collision limit##
###Combust. Flame. 186 (2017) 208-210###

import cantera as ct
import numpy as np
import math

import sys
f = open("Output.out", 'w')
sys.stdout = f

#Mechanism input
ct.add_directory('#####input directory of mechanism########')
gas = ct.Solution('#######Mech_file.cti####################')

#define variables
r = gas.n_reactions
s = gas.n_species
######Change following Values Accordingly#########
T_min = 300    
T_max = 2500
T_step = 50
T = T_min
#v is no. of reactions violating limit and b is no. of bimolecular reactions in mechanism#
v = 0
b = 0

#constants
pi = math.pi
Na = 6.02214129*(10**26)
kb = 1.3806488*(10**(-23))

##########################################
 
##Vioaltion check for forward reactions##

print ('elementary and forward')
print('rxn no.,','rxn,', 'Temperatures')
k_forward = gas.forward_rate_constants
k_eq = gas.equilibrium_constants
stoic_coeff = gas.reactant_stoich_coeff
t = np.zeros(int((T_max - T_min)/T_step)+1)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
 z = i
 for l in range (0,s): tot_coeff = tot_coeff + stoic_coeff(l,i)
 if gas.reaction_type(i)==1 and tot_coeff == 2:
  for j in range(0,s):
   if stoic_coeff(j,i) == 1:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = red_mass_num * gas.molecular_weights[j]
    red_mass_den = red_mass_den + gas.molecular_weights[j]
    coll_dia = coll_dia + diameter
    well_dep = well_dep * well_depth
   elif stoic_coeff(j,i) == 2:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = gas.molecular_weights[j] * gas.molecular_weights[j]
    red_mass_den = 2 * gas.molecular_weights[j]
    coll_dia = 2 * diameter
    well_dep = well_depth * well_depth
  for T in range (T_min, T_max+1, T_step):
   Ts =  kb * T / (math.sqrt(well_dep))
   Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
   k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
   gas.TP = T, 100000
   k_forward[i] = gas.forward_rate_constants[i]
   if k_forward[i] > k_limit:
    if z == i:
     v = v + 1
     print(end='\n')
     print(i, gas.reaction_equations([i]), T , end=' ')
     z = -1
    elif z == -1:
     print (T, end=' ')

  b = b + 1

print(end='\n')
print('# of rxns,', '# of violations')
print (b   , v   )
print (' ')

##########################################

##Vioaltion check for reverse reactions##

print ('elementary and reverse')
print('rxn no.,','rxn,', 'Temperatures')
k_reverse = gas.reverse_rate_constants
stoic_coeff = gas.product_stoich_coeff
v = 0
b = 0
ratiolow=np.zeros(gas.n_reactions)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
 z = i
 for l in range (0,s): tot_coeff = tot_coeff + stoic_coeff(l,i)
 if gas.reaction_type(i)==1 and tot_coeff == 2 and gas.is_reversible(i)==True:
  for j in range(0,s):
   if stoic_coeff(j,i) == 1:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = red_mass_num * gas.molecular_weights[j]
    red_mass_den = red_mass_den + gas.molecular_weights[j]
    coll_dia = coll_dia + diameter
    well_dep = well_dep * well_depth
   elif stoic_coeff(j,i) == 2:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = gas.molecular_weights[j] * gas.molecular_weights[j]
    red_mass_den = 2 * gas.molecular_weights[j]
    coll_dia = 2 * diameter
    well_dep = well_depth * well_depth
  for T in range (T_min, T_max+1, T_step):
   Ts =  kb * T / (math.sqrt(well_dep))
   Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
   k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
   gas.TP = T, 100000
   k_reverse[i] = gas.reverse_rate_constants[i]
   if k_reverse[i] > k_limit:
    if z == i:
     v = v + 1
     print(end='\n')
     print(i,gas.reaction_equations([i]), T, end=' ' )
     z = -1
    elif z == -1:
     print (T, end=' ')

  b = b + 1

print(end='\n')
print('# of rxns,', '# of violations')
print (b   , v  )
print (' ')

#################################################


#falloff rxns#
##########################################

##Vioaltion check for forward reactions##

print ('falloff and forward')
print('rxn no.,','rxn,', 'Temperatures')
stoic_coeff = gas.reactant_stoich_coeff
v = 0
b = 0
ratiolow=np.zeros(gas.n_reactions)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
 z = i
 for l in range (0,s): tot_coeff = tot_coeff + stoic_coeff(l,i)
 if gas.reaction_type(i)==4 and tot_coeff == 2:
  for j in range(0,s):
   if stoic_coeff(j,i) == 1:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = red_mass_num * gas.molecular_weights[j]
    red_mass_den = red_mass_den + gas.molecular_weights[j]
    coll_dia = coll_dia + diameter
    well_dep = well_dep * well_depth
   elif stoic_coeff(j,i) == 2:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = gas.molecular_weights[j] * gas.molecular_weights[j]
    red_mass_den = 2 * gas.molecular_weights[j]
    coll_dia = 2 * diameter
    well_dep = well_depth * well_depth
  for T in range (T_min, T_max+1, T_step):
   Ts =  kb * T / (math.sqrt(well_dep))
   Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
   k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
   ArrA = gas.reaction(i).high_rate.pre_exponential_factor
   EA = gas.reaction(i).high_rate.activation_energy
   Arrb = gas.reaction(i).high_rate.temperature_exponent
   k_forward[i] = ArrA * np.power(T, Arrb) * np.exp(-EA / (T * ct.gas_constant))
   if k_forward[i] > k_limit:
    if z == i:
     v = v + 1
     print(end='\n')
     print(i, gas.reaction_equations([i]), T , end=' ')
     z = -1
    elif z == -1:
     print (T, end=' ')

  b = b + 1

print(end='\n')
print('# of rxns,', '# of violations')
print (b   , v   )
print (' ')

##########################################

##Vioaltion check for reverse reactions##

print ('falloff and reverse')
print('rxn no.,','rxn,', 'Temperatures')
stoic_coeff = gas.product_stoich_coeff
v = 0
b = 0
ratiolow=np.zeros(gas.n_reactions)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
 z = i
 for l in range (0,s): tot_coeff = tot_coeff + stoic_coeff(l,i)
 if gas.reaction_type(i)==4 and tot_coeff == 2 and gas.is_reversible(i)==True:
  for j in range(0,s):
   if stoic_coeff(j,i) == 1:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = red_mass_num * gas.molecular_weights[j]
    red_mass_den = red_mass_den + gas.molecular_weights[j]
    coll_dia = coll_dia + diameter
    well_dep = well_dep * well_depth
   elif stoic_coeff(j,i) == 2:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = gas.molecular_weights[j] * gas.molecular_weights[j]
    red_mass_den = 2 * gas.molecular_weights[j]
    coll_dia = 2 * diameter
    well_dep = well_depth * well_depth
  for T in range (T_min, T_max+1, T_step):
   Ts =  kb * T / (math.sqrt(well_dep))
   Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
   k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
   ArrA = gas.reaction(i).high_rate.pre_exponential_factor
   EA = gas.reaction(i).high_rate.activation_energy
   Arrb = gas.reaction(i).high_rate.temperature_exponent
   gas.TP = T, 100000
   k_eq[i] = gas.equilibrium_constants[i]
   k_forward[i] = ArrA * np.power(T, Arrb) * np.exp(-EA / (T * ct.gas_constant))
   k_reverse[i] = (ArrA * np.power(T, Arrb) * np.exp(-EA / (T * ct.gas_constant))) / k_eq[i]
   if k_reverse[i] > k_limit:
    if z == i:
     v = v + 1
     print(end='\n')
     print(i, gas.reaction_equations([i]), T, end=' ' )
     z = -1
    elif z == -1:
     print (T, end=' ')

  b = b + 1

print(end='\n')
print('# of rxns,', '# of violations')
print (b   , v   )
print (' ')


#################################################


#plog#
##########################################

##Vioaltion check for forward reactions##

print ('plog and forward')
print('rxn no.,','rxn,', 'Pressure,','Temperatures')
stoic_coeff = gas.reactant_stoich_coeff
v = 0
b = 0
ratiolow=np.zeros(gas.n_reactions)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
 for l in range (0,s): tot_coeff = tot_coeff + stoic_coeff(l,i)
 if gas.reaction_type(i)==5 and tot_coeff == 2:
  for j in range(0,s):
   if stoic_coeff(j,i) == 1:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = red_mass_num * gas.molecular_weights[j]
    red_mass_den = red_mass_den + gas.molecular_weights[j]
    coll_dia = coll_dia + diameter
    well_dep = well_dep * well_depth
   elif stoic_coeff(j,i) == 2:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = gas.molecular_weights[j] * gas.molecular_weights[j]
    red_mass_den = 2 * gas.molecular_weights[j]
    coll_dia = 2 * diameter
    well_dep = well_depth * well_depth
  myarray = np.asarray(gas.reaction(i).rates)
  for n in range(0,len(myarray)):
   P = myarray[n,0]
   z = i
   for T in range (T_min, T_max+1, T_step):
    Ts =  kb * T / (math.sqrt(well_dep))
    Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
    k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
    gas.TP = T, P
    k_forward[i] = gas.forward_rate_constants[i]
    if k_forward[i] > k_limit:
     if z == i:
      v = v + 1
      print(end='\n')
      print(i, gas.reaction_equations([i]),P, T, end=' ')
      z = -1
     elif z == -1:
      print (T, end=' ')

  b = b + 1

print(end='\n')
print('# of rxns,', '# of violations(counting violation at each pressure)')
print (b   , v   )
print (' ')

##########################################

##Vioaltion check for reverse reactions##

print ('plog and reverse')
print('rxn no.,','rxn,','Pressure,', 'Temperatures' )
stoic_coeff = gas.product_stoich_coeff
v = 0
b = 0
ratiolow=np.zeros(gas.n_reactions)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
 for l in range (0,s): tot_coeff = tot_coeff + stoic_coeff(l,i)
 if gas.reaction_type(i)==5 and tot_coeff == 2 and gas.is_reversible(i)==True:
  for j in range(0,s):
   if stoic_coeff(j,i) == 1:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = red_mass_num * gas.molecular_weights[j]
    red_mass_den = red_mass_den + gas.molecular_weights[j]
    coll_dia = coll_dia + diameter
    well_dep = well_dep * well_depth
   elif stoic_coeff(j,i) == 2:
    S = gas.species(j)
    diameter = S.transport.diameter
    well_depth = S.transport.well_depth
    red_mass_num = gas.molecular_weights[j] * gas.molecular_weights[j]
    red_mass_den = 2 * gas.molecular_weights[j]
    coll_dia = 2 * diameter
    well_dep = well_depth * well_depth
  myarray = np.asarray(gas.reaction(i).rates)
  for n in range(0,len(myarray)):
   P = myarray[n,0]
   z = i
   for T in range (T_min, T_max+1, T_step):
    Ts =  kb * T / (math.sqrt(well_dep))
    Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
    k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
    gas.TP = T, P
    k_reverse[i] = gas.reverse_rate_constants[i]
    if k_reverse[i] > k_limit:
     if z == i:
      v = v + 1
      print(end='\n')
      print(i, gas.reaction_equations([i]),P, T, end=' ' )
      z = -1
     elif z == -1:
      print (T, end=' ')


  b = b + 1

print(end='\n')
print('# of rxns,', '# of violations(counting violation at each pressure)')
print (b   , v   )
print (' ')

#################################################

#for T in range(T_min, T_max + 1, T_step):

#use this T with the other code

#################################################


f.close()
