#! /usr/bin/python

import sys,math

print 'Deflection'
print 'Evaluation of supersonic flow properties upon deflection'
print 'Daniel Banuti'
print ''


gamma = 1.4


def normalshockrel(M1):
## normal shock relations
  MM1 = M1*M1
  M2 = math.sqrt( (1.0+((gamma-1.0)/2.0)*MM1)/(gamma*MM1-(gamma-1.0)/2.0) )
  p_ratio  = (1.0+(2.0*gamma)/(gamma+1.0)*(MM1-1.))
  T_ratio  = p_ratio*((2.+(gamma-1.)*MM1)/((gamma+1.)*MM1))
  p0_ratio = p_ratio*T_ratio**(gamma/(1.-gamma))
  return M2,p_ratio,p0_ratio,T_ratio


def getbeta(M1,theta):
  MM1 = M1*M1
  delta = 1. #weak shock, set to 0 for strong shock

  Lrad = ((MM1-1.)**2.-3.*(1.+(gamma-1.)/2.*MM1)*((1.+(gamma+1.)/2.*MM1))*(math.tan(theta))**2.)
  if (Lrad > 0.):
    L = math.sqrt( Lrad )
  else:
    print "Angle theta too large for Mach number (lambda). Shock is detached."
    print "Calculate normal shock properties behind detached shock instead."
    print "Lrad="
    print Lrad
    beta = 0.5*math.pi
    return beta
#  print L
  X = ((MM1-1.)**3.-9.*(1.+(gamma-1.)/2.*MM1)*(1.+(gamma-1.)/2.*MM1+(gamma+1.)/4.*MM1*MM1)*(math.tan(theta))**2.)/(L*L*L)
#  print "X=%(X)f" %vars()
  if (abs(X) < 1.):
    beta = math.atan( (MM1-1.+2.*L*math.cos((4.*math.pi*delta+math.acos(X))/3.))/(3.*(1.+(gamma-1.)/2.*MM1)*math.tan(theta) ))
  else:
    print "Angle theta too large for Mach number (Xi). Shock is detached."
    print "Calculate normal shock properties behind detached shock instead."
    print "X Zaehler="
    print ((MM1-1.)**3.-9.*(1.+(gamma-1.)/2.*MM1)*(1.+(gamma-1.)/2.*MM1+(gamma+1.)/4.*MM1*MM1)*(math.tan(theta))**2.)
    print "L^3="
    print (L*L*L)
    print "L="
    print L
    beta = 0.5*math.pi
    return beta
#  print beta
  return beta

def obliqueshockrel(M1,theta,beta=None):
  if beta==None:
  ## if beta not given, it has to be calculated
    beta=getbeta(M1,theta)
    betadeg=math.degrees(beta)
    #print "Beta calculated to %(betadeg)f" %vars()

  Mn1 = M1*math.sin(beta)
  #print "Mn1 calculated to %(Mn1)f" %vars()
  Mn2,p_ratio,p0_ratio,T_ratio = normalshockrel(Mn1)
  M2 = Mn2/math.sin(beta-theta)
  return M2,p_ratio,p0_ratio,T_ratio



def getisentropic(M1,M2):
  MM1 = M1*M1
  MM2 = M2*M2
  T_ratio  = ((1.+(gamma-1.)/2.*MM1)/(1.+(gamma-1.)/2.*MM2)) #T2/T1
  p_ratio  = (T_ratio)**(gamma/(gamma-1.)) #p2/p1
  p0_ratio = 1. #p02/p01
  return p_ratio, p0_ratio, T_ratio

def getnu(M1):
  MM1 = M1*M1
  nu = math.sqrt((gamma+1.)/(gamma-1.))*math.atan(math.sqrt((gamma-1.)/(gamma+1.)*(MM1-1.)))-math.atan(math.sqrt(MM1-1.))
  return nu

def getMfromnu(nu):
  M = iterateM(nu)
  return M

def iterateM(nu_target):
#initial conditions M=2, M=5
  M_1  = 1.2
  M_2  = 2.
  #nu_1 = 0.460413682083-nu_target
  #nu_2 = 1.34251102197-nu_target
  nu_1 = getnu(M_1)-nu_target
  nu_2 = getnu(M_2)-nu_target
#secant loop
  i     = 0
  delta = 1.
  while ((delta > 1.0e-8) and (i < 20)):
    i+=1

    M_0 = M_1
    M_1 = M_2
    nu_0 = nu_1
    nu_1 = nu_2
    #M_2 = M_1 - nu_1*((M_1-M_0)/max(1.0e-8,(nu_1-nu_0)))
    M_2 = M_1 - nu_1*((M_1-M_0)/(nu_1-nu_0))
    #print "M2 = %(M_2)f" %vars()
    nu_2 = getnu(M_2)-nu_target

    delta = abs(nu_2)
    #print "i=%(i)i: delta = %(delta)f" %vars()
  print M_2
  return M_2

def expansion(M1,theta):
  nuM1 = getnu(M1)
  nuM2 = theta + nuM1

  M2 = getMfromnu(nuM2)
  p_ratio, p0_ratio, T_ratio = getisentropic(M1,M2)

  return M2,p_ratio,p0_ratio,T_ratio


def deflect_one(M1,theta,beta=None):
    thetadegs = math.degrees(theta)
    if (abs(theta)< 1.0e-8):
      print "Theta = %(thetadegs)f ~ 0. Evaluating normal shock relations." %vars()
      M2,p_ratio,p0_ratio,T_ratio = normalshockrel(M1)
    else:
      if (theta > 0.):
        print "Theta = %(thetadegs)f > 0. Evaluating oblique shock relations." %vars()
        M2,p_ratio,p0_ratio,T_ratio = obliqueshockrel(M1,theta,beta)
      else:
        print "Theta = %(thetadegs)f < 0. Evaluating Prandtl-Meyer expansion relations." %vars()
        M2,p_ratio,p0_ratio,T_ratio = expansion(M1,abs(theta))

    print "M1 = %(M1)f" %vars()
    print 20*"_"
    print "p2/p1   = %(p_ratio)f" %vars()
    print "T2/T1   = %(T_ratio)f" %vars()
    print "p02/p01 = %(p0_ratio)f" %vars()
    print "M2      = %(M2)f" %vars()
    print
    return M2,p_ratio,p0_ratio,T_ratio


def deflect(M1,theta,beta=None):
  thetadegs = math.degrees(theta)
  if (theta > 0.):
#    print "Theta = %(thetadegs)f > 0. Evaluating oblique shock relations." %vars()
    M2,p_ratio,p0_ratio,T_ratio = obliqueshockrel(M1,theta,beta)
  else:
#    print "Theta = %(thetadegs)f < 0. Evaluating Prandtl-Meyer expansion relations." %vars()
    M2,p_ratio,p0_ratio,T_ratio = expansion(M1,abs(theta))
  return M2,p_ratio,p0_ratio,T_ratio

## MAIN
def computetable():
  outfilename = 'deflection.dat'
  outfile=open(outfilename,mode='w')
  outfile.write('VARIABLES = "M1", "theta", "M2", "p2/p1", "p02/p01", "T2/T1" \n')

  for M1 in range(2,30,1):
    zoneline = 'ZONE T="Mach number = %(M1)f" \n' %vars()
    outfile.write(zoneline)
    for thetadegs in range(-50,50,1):
      theta = math.radians(float(thetadegs))
 #   if (abs(theta)< 1.0e-8):
  #    pass
      #print "Theta = %(thetadegs)f ~ 0. Evaluating normal shock relations." %vars()
      #M2,p_ratio,p0_ratio,T_ratio = normalshockrel(M1)
   # else:
      if (theta > 0.):
        print "Theta = %(thetadegs)f > 0. Evaluating oblique shock relations." %vars()
        M2,p_ratio,p0_ratio,T_ratio = obliqueshockrel(M1,theta,beta)
      else:
        print "Theta = %(thetadegs)f < 0. Evaluating Prandtl-Meyer expansion relations." %vars()
        M2,p_ratio,p0_ratio,T_ratio = expansion(M1,abs(theta))
      outline = str(M1)+'\t'+str(thetadegs)+'\t'+str(M2)+'\t'+str(p_ratio)+'\t'+str(p0_ratio)+'\t'+str(T_ratio)+'\n'
      outfile.write(outline)
  outfile.close()

# M2,p_ratio,p0_ratio,T_ratio = normalshockrel(M1)
#test=obliqueshockrel(M1,p1,theta,2.0)
#test=obliqueshockrel(M1,p1,theta,beta)
#print test

if (len(sys.argv)>1):
  M1        = float(sys.argv[1])
  thetadegs = float(sys.argv[2])
  theta     = math.radians(thetadegs)
  deflect_one(M1,theta)







#def normalshock(