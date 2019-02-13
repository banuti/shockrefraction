#!/usr/bin/python
from  math import *

print'BSD Conditions'
print'author: banu_da (daniel.banuti@dlr.de)'
print'last update: 05.03.2009'
print''



def get_T_from_p_rho(p,rho):
  T=p/(R*rho)
  return T

def get_u_from_M_T(M,T):
  u=M*sqrt(gamma*R*T)
  return u

def get_rho_from_T_p(T,p):
  rho=p/(R*T)
  return rho


################### Gas Properties ##############
# Gas Constant
R = 287.0
# Ratio of Specific Heats
gamma = 1.40



############### Initial Properties #############
# Pressure
p_1=1.0e5 #Pa

# Mach number
M_1=3.0

# Choose arbitrary density
rho_1=1.0 #kg/m^3

# Calculate ideal gas temperature
T_1 = get_T_from_p_rho(p_1,rho_1) #K

# Calculate velocity from speed of sound
u_1 = get_u_from_M_T(M_1,T_1) #m/s




######## Diffraction Layer Properties ############

# Set desired Mach number
M_2 = 1.5

# Equality of velocity and pressure
u_2=u_1
p_2=p_1


# Temperature from Mach numbers
T_2=(M_1/M_2)*(M_1/M_2)*T_1

# Density from ideal gas
rho_2 = get_rho_from_T_p(T_2,p_2)



print " \t 1 \t\t| 2"
print 40*'_'
print "p:   \t %e \t| %e " %(p_1,p_2)
print "M:   \t %e \t| %e " %(M_1,M_2)
print "T:   \t %e \t| %e " %(T_1,T_2)
print "rho: \t %e \t| %e " %(rho_1,rho_2)
print "u:   \t %e \t| %e " %(u_1,u_2)
# Output




print



