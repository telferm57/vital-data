# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 12:44:10 2017

@author: telferm
"""

#%% main pipeline 
pressurefn = 'press.csv.gz'
pressure = getPressureData(pressurefn,None,0)
pressure.head(50)
pressure.describe()
# get the data 
samplerate = 50
smm = getTestData('acc1.csv.gz',nrows=None)
# get onbody  periods within data 
onbodywindow = 2 # minutes window size 
onbodyperiods, offBodyPeriods = dissectData(smm,onbodywindow)

    

# for each period, do the breakdown ... we only have one at the moment,
# so I'll do it manually 
#window size for onbody periods 
for period in onbodyperiods.values():
    startperiod = period[0]*60*onbodywindow*samplerate
    endperiod = (period[1]+1)*60*onbodywindow*samplerate

#pstart = int(startperiod/2)
#pend = int(endperiod/2)
#periodpressure= pressure[0][pstart:pend]
#periodpressure.describe()
#get the vertical acceleration, gravity, body and smoothed versions 
#plt.plot(pressure[0])  

import vitalsensor as vs
gar = vs.getGravity(smm,offBodyPeriods,onbodywindow)

gravity = np.mean(gar)

graverr = (gravity-9.81)/9.81
if graverr > 0.03: print('warning - gravity error > 3%')

ar, arb, argr, ars, vacc, Vhacc = getVectors(smm,startperiod,endperiod,gvalue=gravity)

# check vhacc for coplanarity 
vhnorm = np.cross(Vhacc[0],Vhacc[1])
vhnorm=np.array([vhnorm])
tt = np.dot(vhnorm,Vhacc.transpose())
if np.mean(tt) > 0.05: print('horizontal vactors not co-planar')
smawindow = 50 
sma = getSMA(arb, winsize=smawindow)

actclasses=getActivityClasses(sma,g=gravity)
# sma has one value per window of 1 second. so  timestamps are pos * 50
smat =  np.arange(0,len(sma),1)*50


#get inactive segments - this is also called by sintheta .  

iasegs = getInactiveSegments(sma,5,1.5)
plt.plot(smat,sma)
# get the walks 
walks = getWalks(vacc,actclasses,smawindow)
# get the PTs 
sintheta = getSintheta(ar,sma)
PTlog = getPT(walks,sintheta,samplerate,angleth=0.36)

PTdetail = getPTdetail(vacc,sintheta,PTlog)

#%% file processing : for each subject in file registry, process the file 

  '''  get matlab variables 
    perform daily processing 
    store results in files 
    return registry of results files '''
    
    basedir = 'C:/Users/telferm/projects/vital/Batch1/'
    varsrequired = ['press','indxDn','acc1','realTime']
    for dirs,files in reg.items():
       # if dirs == 'vah002':break
       
       for file in files:
           file = 'vah001_h1_10_00A096236AE0.mat'
           matdata = spio.loadmat(basedir+dirs+'/' + file)
           md = matdata['myData']
           print(file)
           def getmatlabvars(md):
          # for var in varsrequired:
             
               ts = np.ndarray.flatten(md['realTime'][0,0])
               press = md['press'][0,0][:,1] 
               indxDn = np.ravel(md['indxDn'][0,0]) -1

               axyz = md['acc1'][0,0]
               df = pd.DataFrame(axyz,columns=['x','y','z'])
                
               # assemble the pressure data .. 
               pressar = np.ones(len(indxDn))
               indxDn[-1]
               
               for i,j in enumerate(indxDn):
                   pressar[i] = press[j]
                  # if i == 5:break
               df['press'] = pressar
               
               def matlabDate(intime):
                   ''' returns python datetime from matlab timestamp ''' 
                   import datetime
                   days1900 = 719529 # days 1/1/1900 -> 1/1/1970
                   psx = (intime - days1900)*3600.0*24.0
                   t = datetime.datetime.utcfromtimestamp(psx)
                   return t, psx
              
             
               # import 
               ts.transpose().shape
               psxar = [x[1] for x in (map(matlabDate,ts))]
               datetimear = [x[0] for x in (map(matlabDate,ts))]
               #TODO - this is slow - should do in one command 
               
               df['timestamp'] = psxar
               return df
           
           dfsensor = getmatlabvars(md)
           onbodywindow = 2 # minutes window size 
           #TODO standardise all times in seconds 
           onbodyperiods, offBodyPeriods = dissectData(df,onbodywindow)

           
           
           break
          
       break

#%% plotting them  
fig1 = plt.figure()
ax = fig1.add_subplot(111)
#plt.plot(Vvert,label='Vert Vel')

plt.plot(vacc,label='Vert Acc')
Vv = getVertVelocity(vacc)
plt.plot(Vv,label='Vert Vel')
addGrid()
sinfilt = medfilt(sintheta,25)
ax.axhline(y=0.35)
plt.plot(sinfilt, label = 'sin filt 25')
plotWalkTimes(ax,walks,samplerate,{'label':'walks'})
vp1.plotPTTimes(ax,PTlog)
plotInact(ax,iasegs)
vp1.plotPTdetail(ax,PTdetail,50,{})
ax.set_ylim(-5,5)
#x=1;y=3;w=1;h=5
#ax.add_patch(patches.Rectangle((x,y),w,h,fc='b'))
#ax = plt.gca()
ax.minorticks_on()
ax.grid(which='both') 
ax.set_xlabel(r"Samples")
plt.legend()
#ax.plot()
plt.title('Subject 1 - Visit 1')
plt.show()
#%% csv's for Pierre
import csv

featuresdf.to_csv('gaitanalysis1.csv')

with open('activiyADL001.csv', 'w',newline='') as f: 
    w = csv.writer(f)
    w.writerow(['Subject','date','time','sma','actclass'])
    for  row in actdf.iterrows():
       for line in zip(row[1].t,row[1].sma,row[1].actclasses):
             w.writerow([row[1].subject,line[0].strftime('%Y-%m-%d'),
                         line[0].strftime('%H:%M:%S'),
                         line[1],line[2]])
    w.writerows(walks.values())
    

with open('walks010.csv', 'w',newline='') as f: 
    w = csv.writer(f)
    w.writerow(['start_sec','end_sec','steps'])
    w.writerows(walks.values())
    
with open('annotatin006.csv', 'w',newline='') as f: 
    w = csv.writer(f)
    w.writerow(['start_sec','end_sec','steps'])
    w.writerows(annotations.values)
    

ptdet = {i:[v[0],v[1][0]/50,v[2][0]/50,v[3][0]/50] for (i,v) in PTdetail.items()}   

with open('PT010.csv', 'w',newline='') as f: 
    w = csv.writer(f)
    w.writerow(['class','sin','peak','valley'])
    w.writerows(ptdet.values())
    
ptsecs = {i:[v[0]/50,v[1]/50,v[2]/50] for (i,v) in PTlog.items()} 

with open('pts010.csv', 'w',newline='') as f: 
    w = csv.writer(f)
    w.writerow(['start_sec','end_sec','duration'])
    w.writerows(ptsecs.values())
