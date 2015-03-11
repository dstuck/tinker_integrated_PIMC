#! /usr/bin/env python
"""\
Use a simple Discrete Variable Representation method to solve
one-dimensional potentials.

A good general introduction to DVR methods is
Light and Carrington, Adv. Chem. Phys. 114, 263 (2000)
"""

import sys
import re
import numpy as np
import os.path
from optparse import OptionParser
from math import sin,pi,exp,sqrt
from pylab import plot,show,axis
from numpy import zeros,transpose,diag,array
from numpy import dot as matmul
from numpy.linalg import eigh

def ParseInput(ArgsIn):
   '''Parse command line options using outparse.'''
   UseMsg='''usage: %prog [options] directory(default .)
   Check if every input has a completed output file.'''

   parser=OptionParser(usage=UseMsg)
   parser.add_option('-v','--verbose',dest='verbose',action="store_true",help='Print out more information')
   parser.add_option('-t','--tprint',dest='tprint',action="store_true",help='Print out thetas for interpolation')
   parser.add_option('-p','--polynomial',dest='polynomial',action="store_true",help='Use 6th order polynomial fit')
   parser.add_option('-P','--plot',dest='plot',action="store_true",help='Plot fits')
   parser.add_option('-q','--qtype',dest='qtype',type='str',help='Specifies the partion function for the sampling to be done (ho,full,avg). Default = avg')
#   parser.add_option('-H','--hfenergy',dest='hfenergy',action="store_true",help='Store HF energy for O2 job')
#   parser.add_option('-w','--writefile',dest='write',action="store_true",help='Write output to file')
#   parser.add_option('-n','--name',dest='failname',type='str',help='Specifies the name of the output file. Default = fail.txt')
   options,args=parser.parse_args(ArgsIn)
   return options,args

options,args=ParseInput(sys.argv)

def interpolateArray(initArray,nFinal,xMin,xMax,fitOrder=6):
   finalArray = zeros((nFinal,len(initArray[0])))
   theta = zeros((fitOrder+1,len(initArray[0])))
   for n in range(len(initArray[0])):
      nInit = len(initArray[:,n])
      lGrid = xMax[0,n]-xMin[0,n]
      xInit = zeros(nInit)
      for i in range(nInit):
         xInit[i] = xMin[0,n] + i*lGrid/(nInit-1)
      xMat = np.ones((nInit,fitOrder+1))
      xPow = xInit
      for i in range(fitOrder):
         xMat[:,i+1] = xPow
         xPow = xPow*xInit
      theta[:,n] = np.linalg.pinv(xMat).dot(initArray[:,n])
      if options.tprint:
         print "Theta\n",theta[:,n]
#      print "X_guess\n",xMat.dot(theta)
      
      xFinal = zeros(nFinal)
      for i in range(nFinal):
         xFinal[i] = xMin[0,n] + i*lGrid/(nFinal-1)
      xMat = np.ones((nFinal,fitOrder+1))
      xPow = xFinal
      for i in range(fitOrder):
         xMat[:,i+1] = xPow
         xPow = xPow*xFinal

      finalArray[:,n] = xMat.dot(theta[:,n])
   return finalArray,theta

def plotPoly(theta,xMax,nGrid):
   nPlot = len(theta[0])
   nFit  = len(theta)
   #print "nPlot =",nPlot
   #print "nFit =",nFit
   #print "nGrid =",nGrid
   for n in range(nPlot):
      xVec = zeros(nGrid)
      for i in range(nGrid):
         xVec[i] = -xMax[0,n] + i*(xMax[0,n]*2.)/(nGrid-1)
   
      xMat = np.ones((nGrid,nFit))
      xPow = xVec
      for i in range(nFit-1):
         xMat[:,i+1] = xPow
         xPow = xPow*xVec
      
      plot(xVec,xMat.dot(theta[:,n]))
      axis([-xMax[0,n],xMax[0,n],0.0,0.1])
      show()

def getV_1DHO(xVec,w,mass):
   V = zeros(len(xVec))
   for i in range(len(xVec)):
      V[i] = (xVec[i]**2*mass*(w**2)/2.0)
   return V

def readV(filename):
   vlist=[]
   vfile = open(filename,'r')
   line = vfile.readline()
   while(line!=""):
      vlist.append(list(float(str) for str in line.split()))
      line = vfile.readline()
   return np.asarray(vlist)

def readData(filename,offset=0.0):
   vfile = open(filename,'r')
   coords = []
   vList  = []
   line = vfile.readline()
   while(line!=""):
      l = re.search("(\S+?)_(\S+?)_\S+\s+(\S+)",line);
      if(not l==None):
         #print line
         #print l.groups()
         coords.append((int(l.groups()[0]),int(l.groups()[1])))
         vList.append(float(l.groups()[2]))
      line = vfile.readline()
   temp = [0,0]
   for pair in coords:
      temp[0] = max(temp[0],pair[0])
      temp[1] = max(temp[1],pair[1])
   vArray = zeros((temp[0]+1,temp[1]+1))
   for coord, v in zip(coords,vList):
      vArray[coord[0],coord[1]] = v - offset
#   minArray = np.amin(vArray,axis=0)
#   vArray = vArray - minArray
#   print vArray
   #sys.exit()
   return vArray
  
def getT(nGrid,lGrid,mass):
   #print "Getting T"
   T = zeros((nGrid,nGrid))
   for i in range(nGrid):
      for j in range(nGrid):
         if i==j:
            T[i,i] = ((nGrid-1.)*(nGrid-2.)/6.+nGrid/2.)*(pi**2)/mass/lGrid**2
            #T[i,i] = ((nGrid-1.)*(nGrid-2.)/6.+nGrid/2.)/4./mass/lGrid**2
         else:
            T[i,j] = (-1.)**(i-j)/mass*(pi/(lGrid*sin(pi*(i-j)/nGrid)))**2
            #T[i,j] = (-1.)**(i-j)/mass*(1.0/(2.0*lGrid*sin(pi*(i-j)/nGrid)))**2
   return T

def testHarmonic(nGrid,xMin,xMax):
   print "Running Harmonic Test"
   lGrid = xMax - xMin
   #xVec = array[range(xMin,xMax,lGrid/nGrid)]
   xVec = zeros(nGrid)
   for i in range(nGrid):
      xVec[i] = xMin + i*lGrid/(nGrid-1)
   omega = 0.0097741
   mass = 23420.7

   H = getT(nGrid,lGrid,mass)

   V = getV_1DHO(xVec,omega,mass)
   #print xVec
   #print V

   for i in range(nGrid):
      H[i,i] += V[i]
   #print H

   E,U = eigh(H)
   #print "Delta E:"
   #for i in range(nGrid-1):
   #   print (E[i+1]-E[i])
   print E
   #print U
   #print "U",U[0:nGrid,0]

   plot(xVec,V)
#   for i in range(4):
#      plot((xVec[0],xVec[-1]),(E[i],E[i]))
   for i in range(10):
      plot(xVec,U[0:nGrid,i]*E[0]+E[i])
   #show()

   return

def readQChem(xMin,xMax):
   print "Running QChem 1D PES"

   V = readV(args[1])
   
   nGrid = len(V)
   lGrid = xMax - xMin
   xVec = zeros(nGrid)
   for i in range(nGrid):
      xVec[i] = xMin + i*lGrid/(nGrid-1)
   #print xVec
   omega = 0.0097741
   mass = 23420.7

   H = getT(nGrid,lGrid,mass)

   #print xVec
   #print V

   for i in range(nGrid):
      H[i,i] += V[i]
   #print H

   E,U = eigh(H)
   #print E
   #print "Delta E:"
   #for i in range(nGrid-1):
   #   print (E[i+1]-E[i])
   #print U
   #print "U",U[0:nGrid,0]

   plot(xVec,V)
#   for i in range(4):
#      plot((xVec[0],xVec[-1]),(E[i],E[i]))
   for i in range(10):
      plot(xVec,U[0:nGrid,i]*E[0]+E[i])
   show()

   return

def polynomialH(nStates):
   if os.path.isfile(args[1]+".fullQM"):
      V_fullQM = readV(args[1]+".fullQM")
   else:
      V_fullQM = readData("data.txt",readV(args[1]+".equibQM"))
      np.savetxt(args[1]+".fullQM", V_fullQM)
      #outFile = open((args[1]+".fullQM"),'w')
      #outFile.write(np.array_str(V_fullQM).replace('[','').replace(']','')+'\n')
   V_hoQM = readV(args[1]+".hoQM")
   V_fullMM = readV(args[1]+".fullMM")
   V_hoMM = readV(args[1]+".hoMM")
   beta = readV(args[1]+".beta")
   massQM = readV(args[1]+".massQM")
   omegaQM = readV(args[1]+".omegaQM")
   massMM = readV(args[1]+".massMM")
   omegaMM = readV(args[1]+".omegaMM")
   xMax = readV(args[1]+".xMax")
   xMin = -xMax
   if os.path.isfile(args[1]+".nGrid"):
      nGrid = readV(args[1]+".nGrid")
   else:
      nGrid = len(V_fullMM)
   if options.verbose:
      print "interpolating V_fullQM with 6th order polynomial"
   V_fullQM = interpolateArray(V_fullQM,nGrid,xMin,xMax,6)[1]
   for i in range(len(V_fullQM[0])):
      if V_fullQM[6,i]<0:
         print "V_fullQM polynomial goes to -infinity!"
         sys.exit()
   if options.verbose:
      print "interpolating V_fullMM with 6th order polynomial"
   if options.plot:
      plotPoly(V_fullQM,xMax*3.0,nGrid)
      sys.exit()
   V_fullMM = interpolateArray(V_fullMM,nGrid,xMin,xMax,6)[1]
   for i in range(len(V_fullMM[0])):
      if V_fullMM[6,i]<0:
         print "V_fullMM polynomial goes to -infinity!"
         sys.exit()

   nModes = len(omegaQM[0])
   
   vList = [(V_fullQM,V_hoQM,"QM",omegaQM,massQM),(V_fullMM,V_hoMM,"MM",omegaMM,massQM)]
   vDict = {'QM': {'V_full':V_fullQM,'V_ho':V_hoQM,'omega':omegaQM,'mass':massQM},'MM': {'V_full':V_fullMM,'V_ho':V_hoMM,'omega':omegaMM,'mass':massMM}}
   for vType in vDict:
      aV_full   = 0.0
      aV_ho     = 0.0

      for n in range(nModes):
         #print "Mode",n

         H = zeros((nStates,nStates))

         p = -vDict[vType]['omega'][0,n]/2
         r = sqrt(1/vDict[vType]['omega'][0,n]/vDict[vType]['mass'][0,n])
         a = vDict[vType]['V_full'][:,n]
         for i in range(nStates):
            m = float(i)
            #print (-p*(2.0*m+1.0)*0.5)
            #print a[0]
            #print 0.5*a[2]*r**2*(2.0*m+1.0)
            #print 0.25*a[4]*r**4*(6.0*m**2+6.0*m+3)
            #print 0.125*a[6]*r**6*(20.0*m**3+30.0*m**2+40.0*m+15)
            H[i,i] = -p*(2.0*m+1.0)*0.5 + a[0] + 0.5*a[2]*r**2*(2.0*m+1.0) + 0.25*a[4]*r**4*(6.0*m**2+6.0*m+3) + 0.125*a[6]*r**6*(20.0*m**3+30.0*m**2+40.0*m+15)
            if (i+1)<nStates:
               m = float(i+1)
               H[i+1,i] = sqrt(2*m)*(0.5*a[1]*r+0.25*3*a[3]*r**3*m+0.125*a[5]*r**5*(10*m**2+5.0))
               H[i,i+1] = H[i+1,i]
            if (i+2)<nStates:
               m = float(i+2)
               H[i+2,i] = sqrt(m*(m-1))*(0.5*p+0.5*a[2]*r**2+0.25*a[4]*r**4*(4*m-2)+15.0/8*a[6]*r**6*(m**2-m+1))
               H[i,i+2] = H[i+2,i]
            if (i+3)<nStates:
               m = float(i+3)
               H[i+3,i] = sqrt(m*(m-1)*(m-2))*(1/sqrt(8)*a[3]*r**3+sqrt(2)/8*a[5]*r**5*(5*m-5))
               H[i,i+3] = H[i+3,i]
            if (i+4)<nStates:
               m = float(i+4)
               H[i+4,i] = sqrt(m*(m-1)*(m-2)*(m-3))*(0.25*a[4]*r**4+0.125*a[6]*r**6*(6*m-9))
               H[i,i+4] = H[i+4,i]
            if (i+5)<nStates:
               m = float(i+5)
               H[i+5,i] = sqrt(m*(m-1)*(m-2)*(m-3)*(m-4))/sqrt(32)*a[5]*r**5
               H[i,i+5] = H[i+5,i]
            if (i+6)<nStates:
               m = float(i+6)
               H[i+6,i] = sqrt(m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5))*a[6]*r**6/8
               H[i,i+6] = H[i+6,i]

         #print "H:\n",H
   
   #  Diagonalize H
         E,U = eigh(H)
         if options.verbose:
            print "E = ",E[0:10]
         aV_full += E[0]
         aV_ho   += 0.5*vDict[vType]['omega'][0,n]

      print "E_full"+vType+" = ",aV_full
      print "E_ho"+vType+" = ",aV_ho

   return

def deltaAI():
# Expect files to have values listed grid pt X mode (row X col)

   print "Running deltaAI"

   if os.path.isfile(args[1]+".fullQM"):
      V_fullQM = readV(args[1]+".fullQM")
   else:
      V_fullQM = readData("data.txt",readV(args[1]+".equibQM"))
   V_hoQM = readV(args[1]+".hoQM")
   V_fullMM = readV(args[1]+".fullMM")
   V_hoMM = readV(args[1]+".hoMM")
   beta = readV(args[1]+".beta")
   massQM = readV(args[1]+".massQM")
   omegaQM = readV(args[1]+".omegaQM")
   massMM = readV(args[1]+".massMM")
   omegaMM = readV(args[1]+".omegaMM")
   xMax = readV(args[1]+".xMax")
   xMin = -xMax
   if os.path.isfile(args[1]+".nGrid"):
      nGrid = readV(args[1]+".nGrid")
   else:
      nGrid = len(V_fullMM)
   if len(V_fullQM) != nGrid:
      if options.verbose:
         print "interpolating V_fullQM with 6th order polynomial"
      V_fullQM = interpolateArray(V_fullQM,nGrid,xMin,xMax,6)[0]
   if len(V_hoQM) != nGrid:
      if options.verbose:
         print "interpolating V_hoQM with 2nd order polynomial"
      V_hoQM = interpolateArray(V_hoQM,nGrid,xMin,xMax,2)[0]
   if len(V_fullMM) != nGrid:
      if options.verbose:
         print "interpolating V_fullMM with 6th order polynomial"
      V_fullMM = interpolateArray(V_fullMM,nGrid,xMin,xMax,6)[0]
   if len(V_hoMM) != nGrid:
      if options.verbose:
         print "interpolating V_hoMM with 2nd order polynomial"
      V_hoMM = interpolateArray(V_hoMM,nGrid,xMin,xMax,2)[0]

   nModes = len(omegaQM[0])
   
#  Loop over QM and MM to calculate deltaF anh
   vList = [(V_fullQM,V_hoQM,"QM",omegaQM,massQM),(V_fullMM,V_hoMM,"MM",omegaMM,massQM)]
   vDict = {'QM': {'V_full':V_fullQM,'V_ho':V_hoQM,'omega':omegaQM,'mass':massQM},'MM': {'V_full':V_fullMM,'V_ho':V_hoMM,'omega':omegaMM,'mass':massMM}}
   for vType in vDict:
      aV_full   = 0.0
      aV_ho     = 0.0
      aE0       = 0.0
#      aV_fullMM   = 0.0
#      aV_hoMM     = 0.0
#      aV_fullQM   = 0.0
#      aV_hoQM     = 0.0
#      aSL         = 0.0
#      aS          = 0.0
#      aL          = 0.0

      for n in range(nModes):

         lGrid = xMax[0,n] - xMin[0,n]
         xVec = zeros(nGrid)
         for i in range(nGrid):
            xVec[i] = xMin[0,n] + i*lGrid/(nGrid-1)

#         plot(xVec,V_hoQM[:,n])
#         plot(xVec,V_fullQM[:,n])
#         plot(xVec,V_hoMM[:,n])
#         plot(xVec,V_fullMM[:,n])
#         show()
#         sys.exit()

   #  Form H
         H = getT(nGrid,lGrid,vDict[vType]['mass'][0,n])

         #V_avg = 0.25*(V_fullQM[:,n]+V_fullMM[:,n]+V_hoQM[:,n]+V_hoMM[:,n])
         for i in range(nGrid):
            #H[i,i] += V_hoMM[i]                       #0,0
            #H[i,i] += V_fullMM[i]                     #1,0
            #H[i,i] += V_hoQM[i]                       #0,1
            #H[i,i] += V_fullQM[i]                     #1,1
            #H[i,i] += (V_fullQM[i]+V_hoQM[i])*0.5     #1/2,1
            #H[i,i] += (V_fullMM[i]+V_hoMM[i])*0.5     #1/2,0
            #H[i,i] += (V_avg[i])                      #1/2,1/2
            if options.qtype=='full':
               H[i,i] += vDict[vType]['V_full'][i,n]
            elif options.qtype=='ho':
               H[i,i] += vDict[vType]['V_ho'][i,n]
            else:
               H[i,i] += (vDict[vType]['V_full'][i,n]+vDict[vType]['V_ho'][i,n])*0.5     #1/2,1

   #  Diagonalize H
         E,U = eigh(H)
         if options.verbose:
            print "E = ",E[0:10]
         aE0 += E[0]

   #  Set this based on temp
         nMax = min(max(5,int(6/(beta*vDict[vType]['omega'][0,n]))),int(nGrid/2))

   #  Cumulate operators over grid
         cV_full = zeros(nGrid)
         cV_ho   = zeros(nGrid)
#         cV_fullQM = zeros(nGrid)
#         cV_hoQM   = zeros(nGrid)
#         cV_fullMM = zeros(nGrid)
#         cV_hoMM   = zeros(nGrid)
         #cSL       = zeros(nGrid)
         #cS       = zeros(nGrid)
         #cL       = zeros(nGrid)
         cQ        = zeros(nGrid)
         Q         = 0.0
         
         #cV_fullQM = V_fullQM[:].T.dot(U[:,0]**2)
         for i in range(nGrid):
            cV_full   += vDict[vType]['V_full'][i,n]*U[i,:]**2
            cV_ho   += vDict[vType]['V_ho'][i,n]*U[i,:]**2
#            cV_fullQM[i] = V_fullQM.T.dot(U[:,i]**2)
#            cV_hoQM   += V_hoQM[i]*U[i,:]**2
#            cV_fullMM += V_fullMM[i]*U[i,:]**2
#            cV_hoMM   += V_hoMM[i]*U[i,:]**2
#            cSL       += 0.25*(V_fullQM[i]-V_hoQM[i]+V_fullMM[i]-V_hoMM[i])*(V_fullQM[i]+V_hoQM[i]-V_fullMM[i]-V_hoMM[i])*U[i,:]**2
            #cS       += 0.5*(V_fullQM[i]-V_hoQM[i]+V_fullMM[i]-V_hoMM[i])*U[i,:]**2
            #cL       += 0.5*(V_fullQM[i]+V_hoQM[i]-V_fullMM[i]-V_hoMM[i])*U[i,:]**2
            cQ[i]     = exp(-beta*E[i])
            Q         += cQ[i]
         #print V_fullQM[:].T.dot(U[:,0]**2)
         #print cV_fullQM
         if options.verbose:
            print "P(i) =",cQ[0:nMax]/Q

   #  Form averages
         aV_full += cV_full.T.dot(cQ)/Q
         aV_ho += cV_ho.T.dot(cQ)/Q
#         aV_fullMM += cV_fullMM.T.dot(cQ)/Q
#         aV_hoMM += cV_hoMM.T.dot(cQ)/Q
#         aV_fullQM += cV_fullQM.T.dot(cQ)/Q
#         aV_hoQM += cV_hoQM.T.dot(cQ)/Q
#         aSL += cSL.T.dot(cQ)/Q
         #aS += cS.T.dot(cQ)/Q
         #aL += cL.T.dot(cQ)/Q
      print "E_0"+vType+"    =",aE0
      print "<full"+vType+"> =",aV_full
      print "<ho"+vType+">   =",aV_ho
#         print "<fullMM>",aV_fullMM
#         print "<hoMM>",aV_hoMM
#         print "<fullQM>",aV_fullQM
#         print "<hoQM>",aV_hoQM
#         print "<SL>",aSL
      #print "<S>",aS
      #print "<L>",aL
      

#      if(nModes==1):
#         plot(xVec,V_hoQM[:,0])
#         plot([xVec[0],xVec[-1]],[E[0],E[0]])
#         for i in range(10):
#            plot(xVec,U[0:nGrid,i]*E[0]+E[i])
#         show()

   return


if __name__ == '__main__': 
   #testHarmonic(20,-0.2,0.2)
   #readQChem(-0.15,0.15)
   if options.polynomial:
      polynomialH(20)
   else:
      deltaAI()
