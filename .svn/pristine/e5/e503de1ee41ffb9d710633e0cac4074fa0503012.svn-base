#!/usr/bin/env python
import glob
import numpy as np
import os
import sys
import time

#User set variables
#directory = "." #Set directory where profile files are
directory = sys.argv[1]

#Output string to a log file
def output(s,logFile):
  print s
  logFile.write(s)
  logFile.flush()

#Try to open a file, return flag indicating if it exists or not
def openLogFile(path):
  path  = os.path.join(directory, path)
  if os.path.isfile(path):
     print("File exists. Loading.")
     dt = time.strptime( time.ctime(os.path.getctime(path)) )
     fileObj = open(path,'a')
     return fileObj,1,dt
  else:
     print("File does not exist. Creating new one.")
     fileObj = open(path,'w')
     return fileObj,0,None


#Open seeing and xy centres file
seeingFile,contSeeing,stampSE  = openLogFile('seeing.csv')
xyCentFile,contXYCent,stampXY = openLogFile('xyCentre.csv')

#Set stamp to whatever most recent one was loaded, or keep as is
if   stampSE!=None and stampXY!=None: lastStamp = stampXY if stampXY > stampSE else stampSE
elif stampSE!=None: lastStamp = stampSE
elif stampXY!=None: lastStamp = stampXY
else: lastStamp = time.strptime("Mon Sep 1 00:00:00 2018")
lastStamp = time.mktime(lastStamp)

#Print column headers in seeing file if new file
if not contSeeing:
  s = "%20s,"%"DateTime"
  for i in range(8):
    s+="%5s%1i"% ("Star",i)
    if i<7: s+=","
  s+="\n"
  output(s,seeingFile)

#Print column headers in XY file if new file
if not contXYCent:
  s = "%20s," % "DateTime"
  for i in range(8):
    s+="%5s%1i," % ("x",i)
    s+="%5s%i"   % ("y",i)
    if i<7: s+=","
  s+="\n"
  output(s,xyCentFile)
 
while 1:
   
  #Load profiles and extract frame numbers
  profiles = np.array( sorted( glob.glob("%s/*prof.txt" % directory ) ) )
  profiles = np.array([os.path.basename(p) for p in profiles])
  frameNos = np.array( [ int(p.split('_')[1]) for p in profiles] )
  dtStamps = np.array( [ time.strptime( time.ctime(os.path.getctime(os.path.join(directory,p))) ) for p in profiles ] )
  times    = np.array( [ time.mktime(x) for x in dtStamps ])

  #Get only most recent frames (greater than lastFrame)
  use       = times>lastStamp
  useProfs  = profiles[use]
  useNos    = frameNos[use]
  useDTs    = dtStamps[use]
  useTimes  = times[use]

  #Skip and slow down a bit if no new files
  if len(useNos)==0:
    time.sleep(0.5)
    continue

  #If new frames found
  else:

    #Update last frame
    lastStamp = useTimes[-1]

    #Run through new profiles
    for j,prof in enumerate(useProfs):

      #Get sigma from file
      xsig,ysig = 0,0
      for i,line in enumerate(open(os.path.join(directory,prof))):
        if i==2: xc,xsig = (float(v) for v in line.split()[2:4])
        elif i==3: yc,ysig = (float(v) for v in line.split()[2:4])
        elif i>3: break
      sigma = np.sqrt( xsig**2 + ysig**2 )

      #Get frame number and star number
      frame = useNos[j]
      starN = useProfs[j].split('/')[-1].split('_')[0]

      #Get timestamp of file creation time
      dtString = time.strftime("%m/%d/%y %H:%M:%S",useDTs[j])
 
      #Build up strings for new row in each file
      sSE = "%20s," % dtString
      sXY = "%20s," % dtString
      for i in range(8):
        if i==int(starN):
          sSE+="%6.2f" % sigma
          sXY+="%6.2f,%6.2f" % (xc,yc)
        else:
          sSE+="%6s" % "NAN"
          sXY+="%6s,%6s" % ("NAN","NAN")
        if i<7:
          sSE+=","
          sXY+=","
      sSE+="\n"
      sXY+="\n"

      #Output string to each file
      output(sSE,seeingFile)
      output(sXY,xyCentFile)

