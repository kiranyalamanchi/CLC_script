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
######Change following temperatre Accordingly#########
gas.TP = 500, 1.013e+8
T = gas.T
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
print('rxn no.,','rxn,', 'k_f,', 'k_coll')
k_forward = gas.forward_rate_constants
k_eq = gas.equilibrium_constants
stoic_coeff = gas.reactant_stoich_coeff
ratiolow=np.zeros(gas.n_reactions)

for i in range(0,r):
 red_mass_num = 1
 red_mass_den = 0
 coll_dia = 0
 well_dep = 1
 tot_coeff = 0
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
  Ts =  kb * T / (math.sqrt(well_dep))
  Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
  k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
  if k_forward[i] > k_limit:
   v = v + 1
   ratiolow[i] = k_forward[i]/k_limit
   print(i, gas.reaction_equations([i]), "{:.2e}".format(k_forward[i]), "{:.2e}".format(k_limit))
  b = b + 1

print('# of rxns,', '# of violations')
print (b   , v )
print (' ')

##########################################

##Vioaltion check for reverse reactions##

print ('elementary and reverse')
print('rxn no.,','rxn,', 'k_r,', 'k_coll,', 'k_f,', 'k_eq')
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
  Ts = kb * T / (math.sqrt(well_dep))
  Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
  k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
  if k_reverse[i] > k_limit:
   v = v + 1
   ratiolow[i] = k_reverse[i]/k_limit
   print(i, gas.reaction_equations([i]), "{:.2e}".format(k_reverse[i]), "{:.2e}".format(k_limit), "{:.2e}".format(k_forward[i]), "{:.2e}".format(k_eq[i]))
  b = b + 1

print('# of rxns,', '# of violations')
print (b   , v )
print (' ')

#################################################


#falloff rxns#
##########################################
 
##Vioaltion check for forward reactions##

print ('falloff and forward')
print('rxn no.,','rxn,', 'k_f,', 'k_coll')
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
  Ts =  kb * T / (math.sqrt(well_dep))
  Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
  k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
  ArrA = gas.reaction(i).high_rate.pre_exponential_factor
  EA = gas.reaction(i).high_rate.activation_energy
  Arrb = gas.reaction(i).high_rate.temperature_exponent
  k_forward[i] = ArrA*np.power(T,Arrb)*np.exp(-EA/(T*ct.gas_constant))
  if k_forward[i] > k_limit:
   v = v + 1
   ratiolow[i] = k_forward[i]/k_limit
   print(i, gas.reaction_equations([i]), "{:.2e}".format(k_forward[i]), "{:.2e}".format(k_limit))
  b = b + 1

print('# of rxns,', '# of violations')
print (b   , v )
print (' ')

##########################################

##Vioaltion check for reverse reactions##

print ('falloff and reverse')
print('rxn no.,','rxn,', 'k_r,', 'k_coll,', 'k_f,', 'k_eq')
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
  Ts = kb * T / (math.sqrt(well_dep))
  Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
  k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
  ArrA = gas.reaction(i).high_rate.pre_exponential_factor
  EA = gas.reaction(i).high_rate.activation_energy
  Arrb = gas.reaction(i).high_rate.temperature_exponent
  k_forward[i] = ArrA*np.power(T,Arrb)*np.exp(-EA/(T*ct.gas_constant))
  k_reverse[i] = (ArrA*np.power(T,Arrb)*np.exp(-EA/(T*ct.gas_constant)))/k_eq[i]
  if k_reverse[i] > k_limit:
   v = v + 1
   ratiolow[i] = k_reverse[i]/k_limit
   print(i, gas.reaction_equations([i]), "{:.2e}".format(k_reverse[i]), "{:.2e}".format(k_limit), "{:.2e}".format(k_forward[i]), "{:.2e}".format(k_eq[i]))
  b = b + 1

print('# of rxns,', '# of violations')
print (b   , v  )
print (' ')


#################################################


#plog#
##########################################
 
##Vioaltion check for forward reactions##

print ('plog and forward')
print('rxn no.,','rxn,', 'LP,', 'HP,', 'k_f@LP,', 'k_f@HP,', 'k_coll')
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
  Ts =  kb * T / (math.sqrt(well_dep))
  Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
  k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
  myarray = np.asarray(gas.reaction(i).rates)
  for n in range(0,len(myarray)):
   P = myarray[n,0]
   gas.TP = T, P
   k_forward[i] = gas.forward_rate_constants[i]
   if k_forward[i] > k_limit:
    v = v + 1
#    ratiolow[i] = k_forward[i]/k_limit
    print(i, gas.reaction_equations([i]), P, "{:.2e}".format(k_forward[i]), "{:.2e}".format(k_limit))
  b = b + 1

print('# of rxns,', '# of violations(counting violation at each pressure)')
print (b   , v   )
print (' ')

##########################################

##Vioaltion check for reverse reactions##

print ('plog and reverse')
print('rxn no.,','rxn,','LP,', 'HP,', 'k_r@LP,', 'k_r@HP,', 'k_coll,', 'k_f,', 'k_eq')
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
  Ts = kb * T / (math.sqrt(well_dep))
  Om = (1.16145 * (Ts**(-0.14874))) + (0.52487 * (math.exp(-0.7732 * Ts))) + (2.16178 * (math.exp(-2.437887*Ts)))
  k_limit = (math.sqrt((8 * pi * kb * T * red_mass_den * Na) / red_mass_num)) * ((coll_dia/2)**2) * Na * Om
  myarray = np.asarray(gas.reaction(i).rates)
  for n in range(0,len(myarray)):
   P = myarray[n,0]
   gas.TP = T, P
   k_reverse[i] = gas.reverse_rate_constants[i]
   k_forward[i] = gas.forward_rate_constants[i]
   k_eq[i] = gas.equilibrium_constants[i]
   if k_reverse[i] > k_limit:
    v = v + 1
#    ratiolow[i] = k_forward[i]/k_limit
    print(i, gas.reaction_equations([i]), P, "{:.2e}".format(k_reverse[i]), "{:.2e}".format(k_limit), "{:.2e}".format(k_forward[i]), "{:.2e}".format(k_eq[i]))
  b = b + 1

print('# of rxns,', '# of violations(counting violation at each pressure)')
print (b   , v   )
print (' ')

#################################################

f.close()
