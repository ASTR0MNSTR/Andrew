#!/usr/bin/python3
#
#  Tomasz Kwiatkowski           tkastr@vesta.astro.amu.edu.pl
#
#  Showlcs v. 1.5
#
#  2019-05-21 TK, ver 1.3
#
#  * Input data not corrected for light-time nor for unit distances
#    from Sun and asteroid to allow easy identification of the outliers
#  * All plots scaled to the same range in X and Y axes to make comparison easy
#  * JD on the X axis truncated to DD.DDDD so that the mouse cursor coordinate 
#    displayed at the bottom right is shown in full
#  * You can now select ranges of numbers for removal; instead of
#    specifying 1,2,3,16 you can enter 1-3,16
#  * You can now select numbers of the lightcurves to be kept, the
#    rest of lightcurves will be removed
#  * You can use a command line switch "-n" if you do not want the
#    lightcurve plots to b displayed
#
#  2019-11-29 TK, ver. 1.4
#
#  * Dates for the lightcurves are now copied directly from the
#    "OBSERVING TIME:" and refer to the midtime of the lightcurve
#   
#  2020-12-22 TK, ver. 1.5
#   
#  * Improved way the lightcurves are displayed in case of 1, 2, 3 or 4
#    lightcurves
#  * Added a "show_xy" function to display, in the terminal, coordinates
#    of the points which have been clicked upon with a mouse; this can
#    be used to remove the outliers from the respective lighcurves
#
from glob import *
from string import *
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from jdcal import gcal2jd, jd2gcal
from functools import wraps


if len(sys.argv)==1:
   Switch=""
if len(sys.argv)==2:
   Switch=sys.argv[1]

i=-1
docheck=0
Filter = []
PhaseAngle = []
CalDate=[]
AspectData=[]
R=[]
Delta=[]
PA=[]
Lambda=[]
Beta=[]
Ampl=0.0
TimeSpan=0.0
        
def show_xy(event):
    #click x-value
    xdata_click = event.xdata
    ydata_click = event.ydata
    with open('input.txt', "r+") as f:
          old = f.read()  # read everything in the file
          f.seek(0)  # rewind
          f.write(old + (f'{xdata_click} {ydata_click}\n')) #!!!

    print(xdata_click,ydata_click)


#--------------------------------------------------
# Function for reading lightcurves from *.atl file
#--------------------------------------------------
c = 173.145 # speed of light in AU/day
zero_time=0.0
corr="F"

def read_atl_lightcurves(filename):
    LightCurves = []
    global AsteroidName, Ampl, TimeSpan

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    i3=-1
    for line in content:
                if line.startswith("OBJECT........:"):
                        AsteroidName=line.split()[1]
                if line.startswith("LT CORRECTED..:"):
                        corr = line.split()[2].strip()
                if line.startswith("OBSERVING TIME:"):
                        CalDate.append(line.split()[3])
                if line.startswith("ASPECT DATA...:"):
                        d_sun = float(line.split()[2])
                        d_earth = float(line.split()[3])
                        PhaseAngle.append(float(line.split()[4]))
                        R.append(line.split()[2])
                        Delta.append(line.split()[3])
                        PA.append(line.split()[4])
                        Lambda.append(line.split()[5])
                        Beta.append(line.split()[6])                        
                if line.startswith("FILTER........:"):
                        Filter.append(line.split()[1])        
                if line.startswith("DATA:"):
                        lc = []
                        MaxTime=0
                        MinTime=3000000
                        MaxMag=-1000
                        MinMag=1000
                if line.startswith("==============="):
                    lc = np.array(lc)                    
                    LightCurves.append(lc)
                    if Ampl<MaxMag-MinMag:
                       Ampl=MaxMag-MinMag
                    if TimeSpan<MaxTime-MinTime:
                       TimeSpan=MaxTime-MinTime
                    
                try:
                        tokens = line.split()
                        JD = float (tokens[0])
                        mag = float(tokens[1])
                        magsig= float(tokens[2])
                        lc.append([JD, mag, magsig])

                        if MaxTime<JD:
                          MaxTime=JD
                        if MinTime>JD:
                          MinTime=JD

                        if MaxMag<mag:
                          MaxMag=mag
                        if MinMag>mag:
                          MinMag=mag

                except:
                    continue
    return LightCurves

#==============================================================================

def remove_lightcurves(to_remove, filename1, filename2):

#
#   one_curve is a list of lines in one lightcurve 
#   (header and DATA block combined)
#   many curves is a list of lightcurves in *.atl file
#
    one_curve=[]
    many_curves=[]

    f1=open(filename1,'r')
    f2=open(filename2,'w')
    content = f1.readlines()

#
#   read in the content of the *.atl file
#
    i=0
    for line in content:
      if line[:4] != "====":
         one_curve.append(line)
      else:
         many_curves.append(one_curve)
         one_curve=[]
      i=i+1

    for i in range(len(many_curves)):
      printed=0
      remove=0

      for j in range(len(to_remove)):
        k=int(to_remove[j])
        #print("    Curve to be removed= ",k)
        if (i+1)==k:
          remove=1                          

      if remove!=1 and printed!=1:
          f2.write(''.join(many_curves[i]))
          f2.write("===============\n")
          printed=1  
    f2.write("END OF FILE\n")
    
    f1.close()
    f2.close()
#==============================================================================
#
# 2. Ask for input file
#

i=0
print("")
print("")
print("Select the ATL file to work with")
atl_files=glob('*.atl')

for f in atl_files:
  print("%2i    %s" % (i+1, atl_files[i]) )
  i=i+1

Ans=int(input("Enter the number of file to open: "))
print("")
InputFileName=atl_files[Ans-1]
DataFile=InputFileName


#------------------------------------------------
# 3. Read data
#------------------------------------------------
LightCurves=[]
LightCurves = read_atl_lightcurves(InputFileName)
NoOfLightcurves=len(LightCurves)

grid_size= (5,5)
Date=[]
i5=0

print("    Details of lightcurves in the ",InputFileName," file") 
print("===============================================================")
print("No      Date    Filter    r        Delta      PA   lambda  beta")
print("===============================================================")
for lc in LightCurves:
  #Date.append(str(jd2gcal(0.0, float(lc[0,0])[0:3]))
  Y,M,D = CalDate[i5].split('-')
  Y1=int(Y)
  M1=int(M)
  D1=int(D)
  Y1=int(Y)
  M1=int(M)
  D1=int(D)
  PA1=float(PA[i5])
  Lambda1=float(Lambda[i5])
  Beta1=float(Beta[i5])
  print("%2i   %s   %s      %s     %s   %4.1f  %5.1f  %4.1f" %\
       (i5+1,CalDate[i5],Filter[i5],R[i5],Delta[i5],PA1,Lambda1,Beta1) )
              
  i5=i5+1
print("===============================================================")

if Switch != "-n":
  #-------------------------------------------------
  # 4.  Plot lightcurves
  #-------------------------------------------------

  # define the figure size and grid layout properties
  # keep the number of columns fixed at 5, adjust figsize
  #
  i5=0

  if NoOfLightcurves==1:
     cols=1
     rows=1
  if NoOfLightcurves==2:
     cols=2
     rows=1
  if NoOfLightcurves==3:
     cols=3
     rows=1
  if NoOfLightcurves==4:
    cols=4
    rows=1
  if NoOfLightcurves>4:
    cols=5
    rows=NoOfLightcurves // cols + 1

  if rows==4:
    figsize = (16, 8)
  if rows==3:
    figsize = (16, 6)
  if rows==2:
    figsize = (16, 4)
  if rows==1:
    figsize = (16, 3)

  if NoOfLightcurves==1:
    figsize = (4, 3)
  if NoOfLightcurves==2:
    figsize = (8, 3)    
  if NoOfLightcurves==3:
    figsize = (12, 3)
 
  gs = gridspec.GridSpec(rows, cols)
  gs.update(hspace=0.4)

  fig1 = plt.figure(num=1, figsize=figsize)
  ax = []

  for lc in LightCurves:

    col=np.remainder(i5,5)
    row=np.floor_divide(i5,5)
#    if NoOfLightcurves==1:
#      col=1
#      row=1
      
    XMin=lc[0,0]-np.trunc(lc[0,0]/100.0)*100.0
    XMax=XMin+TimeSpan
    AveMag=(np.amax(lc[:,1])+np.amin(lc[:,1]))/2.0
    YMax=AveMag+0.6*Ampl
    YMin=AveMag-0.6*Ampl
    ax.append(fig1.add_subplot(gs[row, col]))
    #ax[-1].get_xaxis().set_visible(False)
    ax[-1].set_xlim([XMin,XMax])
    ax[-1].set_ylim([YMin,YMax])
    #ax[-1].get_yaxis().set_fontsize(20)
    plt.gca().invert_yaxis()
    plt.tick_params(labelsize=8)
    ax[-1].set_title("%3i    %s %s %3.1f%s %s" % (i5+1, CalDate[i5], "PA=", float(PA[i5]), ",", Filter[i5] ), fontsize=8)
    ax[-1].plot( lc[:,0]-np.trunc(lc[:,0]/100.0)*100.0 ,lc[:,1], 'ro')
    i5=i5+1

  plt.suptitle("Lightcurves for %s" % (AsteroidName), fontsize=10)
  plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.9,
                  wspace=None, hspace=None)

  AsteroidPlots=InputFileName.split(".")[0]
  plt.savefig(AsteroidPlots+"_rev.pdf")

  print("Close the graphic window to proceed.")
  print("Its content has been saved in "+AsteroidPlots+"_rev.pdf file")
  
  fig1.canvas.mpl_connect('button_press_event',show_xy)
  plt.show()

#-------------------------------------------------
# 5.  Remove selected LC's
#-------------------------------------------------


print("--------------------------------------------------------------")
print("1) Remove selected lightcurves from ATL file")
print("2) Keep selected lightcurves, remove the rest")
print("3) Quit")
print("--------------------------------------------------------------")
Ans=input("Enter your choice (1 or 2 or 3): ")
if Ans=="1":
  print("Enter the numbers of lightcurves to be removed, separating them")
  print("with a colon. Example: 1,3-8,15")
  Ans=input("Enter the numbers: ")
  Item=Ans.split(",")
  to_remove=[]
  print("Start splitting input string")
  for a in Item:
    print("    Item= ",a)
    if len(a)>2:
       print("    Longer then 2, a range: ",a)
       s=a.split("-")
       for j in range(int(s[0]),int(s[1])+1):
         to_remove.append(str(j))
    else:
       print("    Only one item: ",a)
       to_remove.append(a)
  print(to_remove)  
  print("Enter the name of file for the selected lightcurves")
  print("Example: truncated1.atl")
  ToRemoveFileName=input("Enter the filename: ")
  print(ToRemoveFileName)
  remove_lightcurves(to_remove, InputFileName, ToRemoveFileName)
  print("A new ATL file with selected lightcurves has been created.")
  quit()
if Ans=="2":  
  print("Enter the numbers of lightcurves to be kept, separating them")
  print("with a colon. Example: 1,3-8,15")
  Ans=input("Enter the numbers: ")
  Item=Ans.split(",")
  to_keep=[]
  to_remove=[]
  #print("Start splitting input string")
  for a in Item:
    #print("    Item= ",a)
    if len(a)>2:
       #print("    Longer then 2, a range: ",a)
       s=a.split("-")
       for j in range(int(s[0]),int(s[1])+1):
         to_keep.append(str(j))
    else:
       #print("    Only one item: ",a)
       to_keep.append(a)
  #print(to_keep)
  for i in range(NoOfLightcurves):
      to_remove.append(str(i+1))
  #print("    Old to remove:",to_remove)    
  for b in to_keep:
    to_remove.remove(b)   
  #print("     New to remove:",to_remove)  
  print("Enter the name of file for the selected lightcurves")
  print("Example: truncated1.atl")
  ToRemoveFileName=input("Enter the filename: ")
  print(ToRemoveFileName)
  remove_lightcurves(to_remove, InputFileName, ToRemoveFileName)
  print("A new ATL file with selected lightcurves has been created.")
  quit()
  
if Ans=="3":
  quit()    
