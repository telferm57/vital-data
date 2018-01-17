# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 05:14:55 2017

@author: telferm

analyse clinical test output 

given a participant id : 
    
get the file 

get on body segments - should only be one 

identify start ?

set start to  zero ? (or set annotation start to tod)

activity

try overlay 

walks

pts 

"""
import vitalanalysis.filehandler as fh
import imp
import pandas as pd
import numpy as np
import vitalanalysis.vitalsensor as vs 
import vitalanalysis.vitalplot1 as vp1

import matplotlib.pyplot as plt
from matplotlib import dates as mdates
import imp
from scipy.signal import medfilt
import vitalanalysis.cfg as cfg

#imp.reload(fh)
basedir = cfg.basedir

#imp.reload(vs)
#imp.reload(vp1)


   
#check sample rate 


def magnitude(ar): # assuming x, y, z array 
    ''' given an array with 3 columns or rows, return 1d array of magnitude 
    along the longest axis '''
    arsh = ar.shape  
    if 3 not in arsh:print('I only deal with 3d vectors - or do I ? ')
    if arsh[0]==3: long = 0 
    else: long = 1 
    if long: ar = ar.T
    ar2 = np.square(ar)
    m2 = ar2[0] + ar2[1] + ar2[2]
    m = np.sqrt(m2)
    return m




def countOneGroups(ar):
    ''' will count groupsof more than 1 '1' in array passed in 
    passes back number of groups 
    test : 
    ar =np.array([1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1])
    '''
    import re
    ar = 1*ar # convert to int if it  i'n't already
    np.sum(ar)
    st = ''.join([str(e) for e in list(ar)])
    onegroups = re.findall('1+',st)
    return len(onegroups)


def getSpikes(mag,thresh=50):
    ''' return start and end of 5 - 10 spikes > mag thresh within 3 seconds
    1. get periods of > 50 
    2. find blocks of 3 secs with > 5 of them 
    3. figure it out '''
    spikes = np.greater(mag,thresh)
    spikel = []
    winsz = 200
    overlap = 150
    for  (start, end) in vs.windows(spikes,winsz,overlap):# window size of 5 secs, moving 50 samples 
         window = spikes[start:end]
         if (countOneGroups(window) > 4):
             if  bool(spikel):
                 print('found at ',start)# shrinks all 00 and 010 to 0 
                 if spikel[-1][1] > start:
                     #this is a continuation of previous interval 
                     spikel[-1][1] = start+winsz # update end of previous interval
                 else: # it's a new one
                     spikel.append([start,start+winsz])
             else: # first one 
                 spikel.append([start,start+winsz])
                 
    # TODO delineate end of window by looking at the last spike 
    # in the last window 
    spikel2=[]
    for aspike in spikel:
        # update start to first spike (1) after start 
        spikestart = aspike[0]
        newstart = spikestart + list(spikes[spikestart:]).index(True) # first 1 after spikestart 
        spikel2.append([newstart,aspike[1]])
        
    return spikel2



def fixupAnnot(annotations):
    annotdict = {'sync sensor':'sync',
                 'sync of sensors':'sync',
                 '6m walk normal speed':'6mw normal',
                 '6m walk fast speed':'6mw fast',
                 'stand on one leg right':'stand on R',
                 'stand on one leg left':'stand on L',
            }
   
    annotations['notes2'] = annotations['notes'].map(lambda x: annotdict.setdefault(x.lower(),x))
      
    annotations['notes2'] = annotations['notes2'].map(lambda x: 'natural gait' \
               if x.lower()[0:12] == 'natural gait' else x)
   
    return annotations

def setAnnotStartTime(data,annotations,spikel):
    ''' sets start time (tod)  on annotations matched to the start of the 
    synch point detected in the sensor data '''
    from datetime import timedelta, datetime
    if len(spikel) >2:
        print('fix spikel! start and end max ! ')
        return -1
    #data = newdf.copy()
    # add time index for dataframe x
    annotations['stime']=annotations['start'].map(lambda x : datetime.strptime(x,'%H:%M:%S.%f'))
    # add a small time  to step count times so they do not overlap 
    delta1sec = timedelta(seconds=1)
    annotations['stime'] = annotations.apply(lambda x: x['stime']+ delta1sec if  x['section'] == 'no. steps' else x.stime,
               axis=1)
    
    annotations['delta'] = annotations['stime'].map(lambda x: timedelta(hours=x.hour,
               minutes=x.minute, seconds=x.second, microseconds=x.microsecond))
    # get start time of first synch from signal data 
    startsynch = data.datetime.iloc[spikel[0][0]]
    # get annotations start of synch - there are usually 2 synchs , we want the earliest
    annotsynch = annotations[annotations.notes2=='sync'].stime.min()
    #change it to a delta 
    synchdelta = timedelta(hours=annotsynch.hour,
               minutes=annotsynch.minute, seconds=annotsynch.second,
               microseconds=annotsynch.microsecond)
    # tod = tod from signal + delta from annotations - delta from synch in annotations
    annotations['synched'] = annotations['delta'].map(lambda x: startsynch+x-synchdelta)
    
    #sort annotations 
    return annotations

def adjustPressure(newdf):
    ndf = newdf.copy()
    pressfilt = medfilt(ndf.press.values,25)
    pressmin = ndf.press.min()
    pressmax = ndf.press.max()
    pressmean = ndf.press.mean()
    pressrange=pressmax-pressmin
    ndf['presstd'] = ndf['press'].map(lambda x: ((x-pressmean)/(pressrange/2) ))
    
    return ndf

#%% main for clinical report 
# 

#subject = 'vah001'
#subject = 'vah002'

#subject = 'vah006'
def runClinicalReport(subject,plotmag=False):
    subject = 'vah006'
    basedir = cfg.basedir
    #datasubir = 'Batch1/'
    # rdf is created by create file report in filehandler ... but it is an internal
    # variable for that function. Takes a fair time to build, so it is pickled 
    # this code assumes the pickled file is there. #TODO check file exists 
    #       imp.reload(vs)
    try:
        isinstance(rdf,pd.DataFrame)
    except NameError:
        
        import pickle
     # FileNotFoundError
        report, rdf =  pickle.load(open(cfg.datadir + 'report_rdf.pickle','rb'))
            
           # report, rdf = fh.createFilereport(basedir,datasubir)
           # get clinical monitor data , id off body segment, calculate gravity and
           # off body thresholds 
          # imp.reload(vs)
    dfsensor = vs.getMonitorData(subject,rdf,basedir,adl=0)     
    dfsensor[['x','y','z']].plot()       
    # the following sets the sensor gravity and sensitivity thresholds for 
    # detecting whether onbody or not 
    clinoffsection = {} # these have been derived by observation - speeds up the
    # onbody calculatons below 
    clinoffsection['vah010'] = [90000,120000]
    clinoffsection['vah001'] = [800000,900000]
    clinoffsection['vah002'] = [900000,1000000]
    clinoffsection['vah006'] = [1500000,1600000]
    clinstart = clinoffsection[subject][0]
    clinend = clinoffsection[subject][1]
    global onbodyThresholds
    gravity, thrsd, thrrange, samplerate  = vs.getMonitorThresholds(
                subject,rdf,basedir,clinstart,clinend,dfsensor=False)

    # allow 5% 
    pickle.dump(vs.onbodyThresholds,open(cfg.datadir + 'monitorthresholds',mode='wb'))
    onbodywindow = 2 # minutes window size 
           #TODO standardise all times in seconds 
           
    onBodyPeriods, offBodyPeriods = vs.dissectData(dfsensor,onbodywindow,
                                                   thrsd=thrsd,thrrange=thrrange,
                                                   samplerate=samplerate)
    odsamplerate = {}
    ondata = {} # sensor data corresponding to on body periods 
    for i,period in enumerate(onBodyPeriods.values()):
        #TODO modify for multiple periods - not found necessary yet in clinical data 
        # as the clinical examination has always been the first onbody section
        # probably easier to tie it up with the time ... is that in the clinical results ? 
        startperiod = int(period[0]*60*onbodywindow*samplerate)
        endperiod = int((period[1]+1)*60*onbodywindow*samplerate)
    
        ondata[i] = dfsensor.iloc[startperiod:endperiod,:]
        odsamplerate[i]  =len(ondata[i])/(ondata[i].timestamp.max()-ondata[i].timestamp.min())
        print('smaple rate segment ',i,' ',odsamplerate)
        
    samplerate = odsamplerate[0] # can vary by up to one second    
    
    gravlist = vs.getGravity(dfsensor,offBodyPeriods,onbodywindow,samplerate)
    #vvv = np.asarray(gravlist)
    #vvv = vvv[~np.isnan(vvv)]
    #plt.plot(vvv)
    gmean = np.nanmean(gravlist)
    print('g calculated from whole file: ', gravity)
    print('g calculated from on body section: ', gmean)
    gravity = gmean.copy()
    # note - following vectors are length 7 less than input due to 
    # adjustment for lag after filtering 
    
    ar, arb, argr, ars, vacc, Vhacc, newdf = vs.getVectors(ondata[0],startperiod,endperiod,gvalue=gravity)
    
    epoch = [newdf.datetime.values[0],newdf.datetime.values[-1]] #TODO these are np format - y?
    mag = vs.magnitude(ar)
    if plotmag:
        fig1 = plt.figure()
        plt.plot(newdf.datetime.values,mag)
        fig1.autofmt_xdate()
        plt.show()
        
    spikel = getSpikes(mag,thresh=25) # match with sync from annotations
    #should only be 2 spikes - sometimes others are mis identified 
    if len(spikel) >2:
        print('fix spikel! start and end max ! ')
    if subject in ['vah006','vah001']:
        del  spikel[1]
    if subject in ['vah002']:
        del  spikel[0:2]
        del  spikel[1:3]
       
    annotations = fh.getAnnotations(subject,basedir) #
    annotations = fixupAnnot(annotations)
    annotations = setAnnotStartTime(newdf,annotations,spikel)
    # we now have annoations with timings aligned with the datetime of the sensor data 
    # get walks - first we need the sma values to then  
    smawindow = int(round(samplerate,0)) # first we need the sma values to then...  
    sma = vs.getSMA(arb, winsize=smawindow)
    #TODO check reasonableness for SMA 
    
    actclasses=vs.getActivityClasses(sma,g=gravity) #... compute the activity classes 
    # sma has one value per window of 1 second, hence activity class is per 1 second interval
    
   # iasegs = vs.getInactiveSegments(sma,5,1.5)
    #sinfilt[51380:51390]
    #plt.plot(smat,sma)
    # get the walks 
    walks = vs.getWalks(vacc,actclasses,smawindow,epoch,samplerate,accelthov=1.2)  
    # get postural transitions 
    sintheta = vs.getSintheta(ar,sma) # angle of device with vertical 
    
    PTlog = vs.getPT(walks,sintheta,samplerate,angleth=0.30)
    
    PTdetail = vs.getPTdetail(vacc,sintheta,PTlog,smph=0.30)
    # pressure - create normalised 
    
    return walks, sintheta. PTlog, PTdetail
    

#%% Gait Analysis for each segment in minibest 
def getNearestWalk(stime,walks):
   # stime = annwalks.iloc[0].synched
    # get vector of starttimes
    stimes  =  [x[4] for x in walks.values()]
    # adjust for time differences betwen annotation and walks# stime = annwalk
    stimesdiff = [abs((date -stime).total_seconds()) for date in stimes]
    
    return stimesdiff.index(min(stimesdiff)) + 1 # walks index starts at 1 
   
#imp.reload(vs)
 # imp.reload(mhealthx.extractors.pyGait)
    #subtract this one 
    # get smallest value 
    walklog,arb,walkno='all',samplerate=50,gravity=9.81,plot_test=False
dfgait,walkerrs = vs.getGaitfeatures(walks,arb,walkno='all',samplerate=samplerate)   
gaitcolumns=['subject','test','walkid','start','end',
             'steps1','avg_step_duration','cadence',
             'steps2','sd_step','sd_stride',
             'step_regularity', 'stride_regularity','symmetry'] 
clinical_visit_gait = pd.DataFrame(columns=gaitcolumns)
    
minibest = ['6mw normal','6mw fast','gait speed normal',
            'walk with head turns','walk with pivot turns (walk)','natural gait']

for test in minibest:
    #get start and end
        #get times from annot
       # test = minibest[0]
        annwalks = annotations[annotations.notes2==test]
#        for each one (can be multiple):
        for annwalk in annwalks.synched:
            sensorwalk = getNearestWalk(annwalk,walks)
            print('nearest walk: ',test,annwalk,sensorwalk)
            print('gait report:', dfgait[dfgait.walkid==sensorwalk].number_of_steps)
            
            dfgaitsr = dfgait[dfgait.walkid==sensorwalk] 
            start_time = walks[sensorwalk][4]
            
            if len(dfgaitsr) == 0: 
                gait_row = [subject,test,0,annwalk,start_time,0,
                        0,0,0,0,0,0,0,0]
            else:
                gait_row = [subject,test,int(dfgaitsr.walkid),annwalk,start_time,
                            dfgaitsr.Method1_stepcount[0],
                            dfgaitsr.avg_step_duration[0],dfgaitsr.cadence[0],
                 dfgaitsr.number_of_steps[0],dfgaitsr.sd_step_durations[0],
                 dfgaitsr.sd_stride_durations[0],
                 dfgaitsr.step_regularity[0],dfgaitsr.stride_regularity[0],
                 dfgaitsr.symmetry[0]]
                
            gait_dict = {x:y for (x,y) in zip(gaitcolumns,gait_row)}
            clinical_visit_gait = clinical_visit_gait.append(gait_dict,ignore_index=True)

clinical_visit_gait.to_csv(subject+'clinical_visit_gait.csv')


#%% plot PTs alone 
fig1 = plt.figure()
ax = fig1.add_subplot(111)
 #plt.plot(Vvert,label='Vert Vel')
#create x from datetime column of newdf created manually at the moment above
x = pd.to_datetime(newdf.datetime.values)
#imp.reload(vp1)
x2 = np.array(list(map(mdates.date2num,x)))
Vv = vs.getVertVelocity(vacc)
ax.plot(x2,Vv,label='Vert Vel')
vp1.plotPTTimes(ax,PTlog,epoch) # need to align with annotations, hence x2 
sinfilt = vs.movingAv(sintheta,20)
#plt.plot(x2,sintheta, label = 'sintheta filtered 11')
ax.plot(x2,sinfilt, label = 'sintheta filtered moving average 20',color='g')
#vp1.plotPTdetail(ax,PTdetail,epoch,50,{})
annotax, annottimes  = vp1.pltAnnotations(ax,annotations)


def setylim(lim,ax):
    ax.set_ylim(-ylim,ylim)
    out = ax.plot()
    return out

setylim(2,ax)
ylim=5

annotax.set_ylim(-ylim,ylim)
ax.set_ylim(-ylim,ylim)
annotax.set_xlim(annottimes[0],annottimes[-1])
ax2.set_xlim(annottimes[0],annottimes[-1]) 
ax.set_xticks(annottimes)
hours = mdates.HourLocator()   # every hour
minutes2 = mdates.MinuteLocator(interval=2)  # every month
seconds = mdates.SecondLocator(interval=5)
hoursFmt = mdates.DateFormatter('%H:%M')
secsFmt = mdates.DateFormatter('%H:%M:%S')
# format the ticks
ax.xaxis.set_major_locator(seconds)
ax.xaxis.set_major_formatter(secsFmt)
ax.xaxis.set_minor_locator(minutes)

annottimes=mdates.date2num(annotations.synched.tolist())
ax.set_xlim(annottimes[0],annottimes[-1])
ax.minorticks_on()
print(type(ax.get_xmajorticklabels())) 
ax.grid(which='both')
annotax.grid(which='both')
    #ax2.set_xlim(tst[0],tst[-1])
dd = annotax.get_yticklabels().label.set_fontsize(6) 
zed = [tick.label.set_fontsize(6) for tick in annotax.xaxis.get_major_ticks()]

timeFmt = mdates.DateFormatter('%M:%S')
ax.xaxis.set_major_formatter(timeFmt)
newdf = adjustPressure(newdf)
ax.plot(x2,newdf.presstd)
#plt.show()
ax.legend()
#labels = ax.get_xticklabels() - does not work ?!
#plt.setp(labels, rotation=30, fontsize=10)
fig1.autofmt_xdate()
# attempt t add another axis for samples rather than time  
ticloc = np.arange(0,len(newdf),3000)
ax2 = fig1.add_axes((0.1,0.1,0.8,0.0))
ax2.yaxis.set_visible(False) # hide the yaxis
ax2.set_xticks(ticloc)
plt.show()
#%% plot walks alone
#imp.reload(vp1)
fig1 = plt.figure()
ax = fig1.add_subplot(111)
#plt.plot(Vvert,label='Vert Vel')
#create x from datetime column of newdf created manually at the moment above
x = pd.to_datetime(newdf.datetime.values)

x2 = np.array(list(map(mdates.date2num,x))) # convert to  matplotlib dates 
Vv = vs.getVertVelocity(vacc)
#ax.plot(x2,Vv,label='Vert Vel')

ax.plot(x2,vacc,label='Vert Acc')
ax.set_title(subject)
vp1.plotWalkTimes(ax,walks,epoch,samplerate,5,{'label':'walks'})
timeFmt = mdates.DateFormatter('%d/%m %H:%M:%S')
ax.xaxis.set_major_formatter(timeFmt)
fig1.autofmt_xdate()
annotax,annottimes  = vp1.pltAnnotations(ax,annotations)
ax.set_xlim(annottimes[0],annottimes[-1])
#newdf = adjustPressure(newdf)
#plt.plot(x2,newdf.presstd)
plt.show()
ax.legend()
#%% plotting walks and pts
import vitalplot1 as vp1
imp.reload(vp1)
fig1 = plt.figure()
ax = fig1.add_subplot(111)
#plt.plot(Vvert,label='Vert Vel')
#create x from datetime column of newdf created manually at the moment above
x = pd.to_datetime(newdf.datetime.values)

x2 = np.array(list(map(mdates.date2num,x)))

plt.plot(x2,vacc,label='Vert Acc')
Vv = vs.getVertVelocity(vacc)
plt.plot(x2,Vv,label='Vert Vel')
vp1.addGrid()
vp1.addGrid()
vp1.plotWalkTimes(ax,walks,epoch,samplerate,10,{'label':'walks'})
timeFmt = mdates.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(timeFmt)
fig1.autofmt_xdate()
vp1.plotPTTimes(ax,PTlog,PTdetail,epoch) # need to align with annotations, hence x2 
vp1.plotInact(ax,iasegs)
#sinfilt = medfilt(sintheta,11)
#for x in [k for k,v in PTdetail.items() if v == -1]:
#    PTdetail.pop(x)
    

sinfilt = vs.movingAv(sintheta,20)
#plt.plot(x2,sintheta, label = 'sintheta filtered 11')
plt.plot(x2,sinfilt, label = 'sintheta filtered moving average 20',color='g')
vp1.plotPTdetail(ax,PTdetail,epoch,50,{})
annotax,annottimes  = vp1.pltAnnotations(ax,annotations)
ax.plot()
plt.plot(x2,newdf.presstd.values,label='pressure') 
plt.show()
plt.legend()

#%% read and plot annotation data 
#TODO put timestamps on these                
smafiles,vectorfiles = multifileSMA(reg,basedir,subjects=['vah002'])

activityLevels = multifileActclasses(smafiles,gravity=9.81)

#times = mdates.date2num(tst.values)
#convert synched  annot times to matplotlib 

fig, ax = plt.subplots()



# plot vacc for trial  
plt.plot(newdf['datetime'].values,vacc)
plt.axhline(y=0)
plt.axhline(y=1.4,c='r')
plt.plot([],[])
#plt.scatter(tst.values,np.random.random(69))

#fig.autofmt_xdate() too complex
ax.minorticks_on()
# dont use grid use vertical lines 
ax.grid(which='both')
#plt.gcf().set_size_inches(32,24)

fig,ax = plt.plot()

        
#%% plot 'monitor active' times  for all participants       

from datetime import time,datetime, date

subjects = cfg.subjects
# ASSUMING we have all the actdf data pickled in the subject directory 
# create the onbody data frame for use in graphic below 
ADLactsummary = pd.DataFrame([])
td2 = timedelta(hours=2)

for subject in subjects:
    # subject = 'vah002'
    ADLactdf = fh.getACTdf(basedir,subdir,subject)
    ADLactdf['date']=ADLactdf.t.map(lambda x: x[0].date()) 
    ADLactdf['dateend']=ADLactdf.t.map(lambda x: x[-1].date())        
    ADLactdf['start']=ADLactdf.t.map(lambda x: x[0].time())        
    ADLactdf['end']=ADLactdf.t.map(lambda x: x[-1].time()) 
    ADLactdf['dayoffset'] = ADLactdf.date.map(lambda x: datetime(x.year,x.month,x.day)+td2)
    # remove sma and activity values from summary - they are no longer meaningful

    ADLactdf.drop(labels=['t','actclasses','sma'],inplace=True,axis=1)
    # if any records span 2 days:
    tt2 = pd.DataFrame([])
    tt = ADLactdf[ADLactdf.date != ADLactdf.dateend]
    ADLactdf = ADLactdf.drop(tt.index)
    print('number of rows spaning > 1 day',subject,len(tt))
    
    for i,rec1 in tt.iterrows():
#        stdate = datetime.strptime(rec1.date,'%Y-%m-%d')
#        edate = datetime.strptime(rec1.dateend,'%Y-%m-%d')
        firstrec = True

        stdate = rec1.date
        edate = rec1.dateend
        dif = edate - stdate
        numrecs = dif.days+1
        finalendtime = rec1.end
        for i in np.arange(numrecs):
            
            rec1.end ='23:59:59'
       
            nsdate = stdate + timedelta(days=int(i))
            #nsdate = datetime.strftime(nsdate,'%Y-%m-%d')
            rec1.date = nsdate # will be 0 for 1st
            rec1.dateend = nsdate
            if not firstrec:
                rec1.start='00:00:00'
            if i == numrecs-1: # last iteration - overwrite end time with end time from original 
                 rec1.end = finalendtime
                 print('last')
            firstrec = False
            tt2 = tt2.append(rec1)
            #print(rec1)    
            
    ADLactdf = pd.concat([ADLactdf,tt2],axis=0,
                        ignore_index=True)
    
    ADLactdf['dayoffset'] = ADLactdf.date.map(lambda x: datetime(x.year,x.month,x.day)+td2)
    ADLactsummary = pd.concat([ADLactsummary,ADLactdf],axis=0)
         
ADLactsummary = ADLactsummary.sort_values(['subject','date','start'])
#plt.plot_date(x,y,ydate=True)
#plt.plot_date(x,end,ydate=True,c='r')
#imp.reload(plt)

# reconfigure rdf to allow for files  > 1 days where 
# iterate through rdf df 
# if startdate != end date, add new records for # days 
#   else add this recrd to new df 
rdf2 = pd.DataFrame([])
for ix,ser in rdf.iterrows():
    if ser.date == ser.dateend:
        rdf2 = rdf2.append(ser)
        type(rdf.get('date'))
    else:
        # get dif in dates 
       # ser = rdf.iloc[21]
        stdate = datetime.strptime(ser.date,'%Y-%m-%d')
        edate = datetime.strptime(ser.dateend,'%Y-%m-%d')
        dif = edate - stdate
        # set existing record to end time 23:59:59 and 
        # enddate same date 
        # create new record with next day start date and
        
        # 23:59:59 end or existing  end time if this is the last 
        # record to add 
        rec1 = ser.copy()
        firstrec = True
        for i in np.arange(dif.days+1):
          #  i=0
            rec1 = ser.copy()
            rec1.end ='23:59:59'
            rec1.dateend = rec1.date
            nsdate = stdate + timedelta(days=int(i))
            nsdate = datetime.strftime(nsdate,'%Y-%m-%d')
            rec1.date = nsdate # will be 0 for 1st
            rec1.dateend = nsdate
            if not firstrec:
                rec1.start='00:00:00'
            if i == dif.days: # last iteration - overwrite end time with ser endtime 
                 rec1.end = ser.end
                 print('last')
            firstrec = False
            rdf2 = rdf2.append(rec1)
            #print(rec1)    
                 
        
rdf2 = rdf2[rdf2.adl==1]  # remove non-adl records 



def plotVert(ax1,dates,starttime,endtime,subject,kwargs):
    print(kwargs)
    print('d')
    for j in range(len(dates)):
        ax1.plot_date([dates[j],dates[j]],[starttime[j],endtime[j]],**kwargs)
        ax1.set_title(subject)
        ax1.xaxis.set_major_formatter(timeFmt)
    return


axn = [0,0,0,0]
f1, ((axn[0], axn[1]), (axn[2], axn[3])) = plt.subplots(2, 2,  sharey='row')
timeFmt = mdates.DateFormatter('%d/%m')

for i,subject in enumerate(subjects):
    #subject='vah001'
    x = rdf2[rdf2.subject==subject].date.apply(lambda x: datetime.strptime(x,'%Y-%m-%d').date()).values
    y = rdf2[rdf2.subject==subject].start.apply(lambda x: datetime.strptime(x,'%H:%M:%S').time()).values
    end = rdf2[rdf2.subject==subject].end.apply(lambda x: datetime.strptime(x,'%H:%M:%S').time()).values
    
    kwargs=dict(linestyle='-',c='r',linewidth=2.0)
    plotVert(axn[i],x,y,end,subject,kwargs)
    x=ADLactsummary[ADLactsummary.subject==subject]['dayoffset'].values
  
    y = ADLactsummary[ADLactsummary.subject==subject]['start'].values
    end = ADLactsummary[ADLactsummary.subject==subject]['end'].values
    kwargs=dict(linestyle='-',c='g',linewidth=1,marker='+')
    plotVert(axn[i],x,y,end,subject,kwargs)
        
# from adlact data frame get on body periods in same format as rdf  - sdate,date end, start, end 
 
x=ADLactsummary['date'].values
y = ADLactsummary['start'].values
end = ADLactsummary['end'].values

plotVert(axn[2],x,y,end,linetype='go-')     
       
yy = ADLactdf.t.values[0]   
f1.suptitle(' monitor wear times')
    
fig1.autofmt_datesdate()

# get one actvity line for this subject .. sma and time from actdf 
# which is generated in vitalADL.manageActvitylevels and contains all sma and rdf
# for on-body protions of files 

tt = actdf.iloc[3]
#now plot as a vertical color map on time axis 

tt.index
ts = tt.t
sma = tt.sma
len(ts) == len(sma)
plt.plot(x=1,y=sma)

plt.axvspan(x[3],x[4],fc='b')
#plot acivity levels 

data = np.arange(100, 0, -1).reshape(10, 10)

fig, ax = plt.subplots()
cax = fig.add_axes([0.27, 0.8, 0.5, 0.05])

im = ax.imshow(data, cmap='gist_earth')
fig.colorbar(im, cax=cax, orientation='horizontal')
plt.show()
a1 = activityLevels['vah002']['vah002h2800A096236B0F.mat']

x = np.arange(6)
y = np.arange(5)
z = x * y[:, np.newaxis]

#%% sample heatmap based on date/time
# generate a range of times over 2 dates 
tt = []
for h in range(8):
    #print(type(tt))
    for m in np.arange(0,46,15):
        tt.append(datetime(2017,1,18,h,m,0))
len(ff)
ff = np.random.random(32) # generate  random values  
dftt = pd.DataFrame(data=ff,index=tt)   
 
# get unique dates (i.e. date/month/year) from data 
dddates = np.array([date.isoformat(x) for x in dftt.index.date])
dddates = list(set(dddates))
dddates.sort(reverse=True)

# put different date's data into columns 
datear={}   
for dtsa in dddates:
   s1 = dftt[dtsa][0] # get series 
   s1.index = s1.index.time
   datear[dtsa] = s1
   

dtff2 = pd.DataFrame.from_dict(datear)
# plot it 
# https://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor
fig, ax = plt.subplots()
cmap = plt.get_cmap('Blues')
# see https://stackoverflow.com/questions/9214971/cmap-set-bad-not-showing-any-effect-with-pcolor
# for handling nans in the data 
cmap.set_bad(color = 'k', alpha = 1.)
heatmap = ax.pcolormesh(dtff2, cmap=cmap, alpha=0.8)
# make some cells nan (where there is no data )
dtff2.loc[(dtff2.index < time(5,30,0)) & (dtff2.index > time(2,30,0) ) ]

# put the major ticks at the middle of each cell

ax.set_yticks(np.arange(dtff2.shape[0]) + 0.5, minor=False)
ax.set_xticks(np.arange(dtff2.shape[1]) + 0.5, minor=False)
ax.set_xticklabels(dtff2.columns, minor=False)
ax.set_yticklabels(dtff2.index, minor=False)


