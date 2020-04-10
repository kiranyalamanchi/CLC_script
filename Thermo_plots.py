import cantera as ct
import numpy as np
import math
import sys

f = open("test.out", 'w')
sys.stdout = f

#Mechanism input
ct.add_directory('\yalamakk\Desktop\Kiran_intern_data\Kiran_ios\CloudFlame')
gas = ct.Solution('Aramco1.3_mech.cti')

#define variables
s = gas.n_species
i = 115

def thermocalculator(stoic_coeff):
   for j in range(0,s):
       if stoic_coeff(j, i) != 0:
           S = gas.species(j)
           T_start = S.thermo.min_temp
           T_max = S.thermo.max_temp
           print (S, T_start, T_max)
           k = open('{0}.txt'.format(j),'w')
           k.write("#This file has data for {}\n".format(S))
           for T in range(int(T_start), int(T_max + 1), 25):
               Cp = S.thermo.cp(T)
               enthalpy = S.thermo.h(T)
               entropy = S.thermo.s(T)
               k.write ("{},{},{},{}\n".format(T, Cp, enthalpy, entropy))

stoic_coeff = gas.reactant_stoich_coeff
thermocalculator(stoic_coeff)

stoic_coeff = gas.product_stoich_coeff
thermocalculator(stoic_coeff)



