#! /usr/bin/python
# -*- coding: utf-8 -*-


import sys,math,deflection

print "Calculate thermally stratified ramp flow"


def getdp(M2,M3,dtheta,theta1,p_31):
  degdtheta = math.degrees(dtheta)
  M4,p_43,p0_43,T_43 = deflection.deflect(M3,-dtheta)
  degdtheta = math.degrees(dtheta)
  theta2 = theta1+dtheta  
  theta2deg = math.degrees(theta2)
  M5,p_52,p0_52,T_52 = deflection.deflect(M2,theta2)
  p_41 = p_31*p_43
  delta_p = p_52-p_41
  return delta_p,p_52,p0_52,p0_43,p_41,M4


def iterate_dtheta(M1,M2,theta1):
#  print "Start iteration."
  M3,p_31,p0_31,T_31 = deflection.deflect(M1,theta1)
  if (M3<1.):
#    print "M3 is subsonic."
    return -6.,-10.,-10.,-10.,-10.,-10.
  #initial conditions
  dtheta_1  = math.radians(0.)
  dp_1,p_52,p0_52,p0_43,p_41,M4 = getdp(M2,M3,dtheta_1,theta1,p_31)
  print "Init dp1 = %(dp_1)f" %vars()
  dtheta_2  = math.radians(5.)
  dp_2,p_52,p0_52,p0_43,p_41,M4 = getdp(M2,M3,dtheta_2,theta1,p_31)
#  print "Init dp2 = %(dp_2)f" %vars()
  #secant loop
  i     = 0
  delta = 1.
  while ((delta > 1.0e-8) and (i < 10)):
    i+=1
    dtheta_0 = dtheta_1
    dtheta_1 = dtheta_2
    dp_0 = dp_1
    dp_1 = dp_2
    degdtheta = math.degrees(dtheta_2)
    dtheta_2 = dtheta_1 - dp_1*((dtheta_1-dtheta_0)/(dp_1-dp_0))
    degdtheta = math.degrees(dtheta_2)
    dp_2,p_52,p0_52,p0_43,p_41,M4 = getdp(M2,M3,dtheta_2,theta1,p_31)
    degdtheta = math.degrees(dtheta_2)
#    print "=======  dtheta = %(degdtheta)f" %vars()
    delta = abs(dp_2-dp_1)
#    print "=======  i = %(i)i: |dp| = %(delta)f" %vars()
  dtheta = dtheta_2
  p0_41 = p0_43*p0_31
  M6,p_64,p0_64,T_64 = deflection.deflect(M4,-dtheta)
  p_61  = p_64*p_41
  p0_61 = p0_64*p0_41
  return dtheta,p_52,p0_52,p0_41,p_61,p0_61


 
def computethetarange(M1,M2,theta1):
  #initial ramp oblique shock
  M3,p_31,p0_31,T_31 = deflection.deflect(M1,theta1)
#  print "p3/p1 = %(p_31)f" %vars()
#  print "M3 = %(M3)f" %vars()

  outfilename = 'refractionM'+str(M1)+'M'+str(M2)+'.dat'
  outfile=open(outfilename,mode='w')
  outfile.write('VARIABLES = "dtheta", "p5/p2", "p4/p1", "p3/p1", "p4/p3", "dp" \n')

  thetastep=0.01
  mindegs = -40
  maxdegs = 40
  minrange = int(mindegs/thetastep)
  maxrange = int(maxdegs/thetastep)
  for dthetadegs_it in range(minrange,maxrange,1):
    dthetadegs=thetastep*dthetadegs_it
    dtheta = math.radians(float(dthetadegs))
  #bottom: shock or expansion
    M4,p_43,p0_43,T_43 = deflection.deflect(M3,-dtheta)
    theta2 = theta1+dtheta
#    print "theta2 %(theta2)f = theta1 %(theta1)f + dtheta %(dtheta)f" %vars()

#    print "p4/p3 = %(p_43)f" %vars()
    M5,p_52,p0_52,T_52 = deflection.deflect(M2,theta2)
#    print "p5/p2 = %(p_52)f" %vars()
    p_41 = p_31*p_43
#    print "p4/p1 = p3/p1 %(p_31)f * p4/p3 %(p_43)f" %vars()
#    print "p3/p1 = %(p_31)f" %vars()
    dp=p_52-p_41
    outline = str(dthetadegs)+'\t'+str(p_52)+'\t'+str(p_41)+ '\t'+str(p_31)+'\t'+str(p_43)+'\t'+str(dp)+'\n'
    outfile.write(outline)
  outfile.close()
  return 


def computethetamap():
  outfilename = 'thetaMach.dat'
  outfile=open(outfilename,mode='w')
  outfile.write('VARIABLES = "M1", "M2", "dtheta", "p5/p2", "p05/p02", "p04/p01", "p6/p1", "p06/p01" \n')

  maxM  = 10.025
  minM  = 1.025
  stepM = 0.025

  minMint = int(minM/stepM)
  maxMint = int(maxM/stepM)

  noM = maxMint-minMint

  for theta1deg in range(10,25,5):
    zoneline = 'ZONE T="Theta 1 = %(theta1deg)f°" I=%(noM)i, J=%(noM)i, F=POINT \n' %vars()
    outfile.write(zoneline)
    print zoneline
    theta1 = math.radians(theta1deg)
    for M1it in range(minMint,maxMint,1):
      for M2it in range(minMint,maxMint,1):
        M1=M1it*stepM
        M2=M2it*stepM
        dtheta,p_52,p0_52,p0_41,p_61,p0_61 = iterate_dtheta(M1,M2,theta1)
        dthetadeg = math.degrees(dtheta)
        outline = str(M1)+'\t'+str(M2)+'\t'+str(dthetadeg)+'\t'+str(p_52)+'\t'+str(p0_52)+'\t'+str(p0_41)+'\t'+str(p_61)+'\t'+str(p0_61)+'\n'
        outfile.write(outline)
  outfile.close()
  return



M1      = 6.
M2      = 10.
theta1  = math.radians(20.)

computethetarange(M1,M2,theta1)

#dtheta,p_52,p0_52,p0_41 = iterate_dtheta(M1,M2,theta1)
#dthetadeg = math.degrees(dtheta)
#print "dtheta = %(dthetadeg)f" %vars()

#computethetamap()





