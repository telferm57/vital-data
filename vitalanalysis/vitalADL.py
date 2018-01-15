# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:31:20 2017

@author: telferm

basic workflow handling for sensor files. Once the basics are identified, 
look at implementing in a python pipline like luigi or airflow 

input to this process is a registry of files that comprise a number of days of
sensor data , perhaps with gaps of time and mutiple files per day . 
we then create a registry of time periods where the device is being worn 
('on body'). These are the time periods that will be analysed

Do we create a new set of files that are only on body activity ? I think not

We could create one file per day though , that would simplify. 
then run 'plugins' to perform each analysis.

plugins can be sensor processors, aggregators or visualisers 

format of sensor plugin: 
    
    input: 
        sensor file with timestamps and sensor data 
        OR output of previous plugin , such as window values for SMA 
        
        
        descriptor of sensor data - must map the following : 
            Gyro axes: gyrox, gyroy, gyroz 
            pressure: Press ,mmHG 
            Accelerometer: x,y,z m/sec , gravity included 
        
        prerequisties: e.g. Gait analysis requires walk segments and activity 
            levels 
            
    output: 
        report file which is suitable for visualiser andor aggregator plugins 
        
        
    

"""

''' sample first plugin - SMA  '''


import vitalsensor as vs
import filehandler as fh
import vitalutils as vu
import imp
import pandas as pd
import pickle 
import numpy as np
#imp.reload(vu)
#imp.reload(vs)
#imp.reload(mhealthx.extractors)
vs.initValues()
import vitalplot1 as vp1
import matplotlib.pyplot as plt
from matplotlib import dates as mdates
from datetime import datetime,timedelta
def pluginSMA(data,datadescriptor,prereq,windowsz):
    ''' anadoning for now 13/11 11:40 - taking too much time away ! ''' 
    
def pluginActivityLevel(data,datadescriptor,prereq): # obv will be class 
    ''' activity level needs only accelerometer data and SMA (signal magnitude
    area) to generate a classification of activity in timeslots 
        data: accelerometer data 
        datadescriptor : this is actually SMA for a bunch of windows that relates 
        to a partiular file 
        
    '''
    # check data matches expected data descriptor
    # check prereq available 
#%% acitivity processing 
def manageActvitylevels(subject,rdf,basedir,subdir,walksgait=False):
    ''' for an individual subject, for each adl file, create an activity level
    file and a registry of those files.
    
    input : 
        subject: id of trial particiant
        rdf : list of file names per subject
        basedir: base location of files 
        subdir : subdirectory of location of adl files 
            (it is assumed each subject has a ditrectory named as the subject id )
            
     output:
         actfilereg: registry of files containing activty data per subject
             location: subdirectory 'adl' of subject in basedir/subdir 
         activity file: subject, timestamp, activity level fine, 
         activity level catagorised 
        pseudo:
            for each file in rdf:
                get on body/inactive segments 
                get activity levels 
                generate file 
                generate dict entry for actfilereg
                '''
    basedir = 'C:/Users/telferm/projects/vital/data/'
    subdir = 'Batch1/'
   
    walksgait=False
    threshsd = 0.015 # standard deviation threshold  of block of 25 steps 
   # del rdf

    adlfiles =  rdf[(rdf.subject==subject) & (rdf.adl==1)][['fn']]
    adlfilestime =  rdf[(rdf.subject==subject) & (rdf.adl==1)][['date','start']]
    #TODO put this logic in rdf 
    adlfilestime['startdatetime'] =                                     \
    adlfilestime.apply(lambda x: datetime.strptime(x['date'] +          \
                                                   'T' + x['start'],    \
                                                   '%Y-%m-%dT%H:%M:%S'),axis=1)
    
    adlfiles = adlfiles.values.ravel()
    activitydata = {}
    samplerate=50 #TODO derive from data 
    ADLactdf = pd.DataFrame({}) # activities and sma stored here
    ADLwalkdf = pd.DataFrame({}) # walk details 
    ADLgaitdf = pd.DataFrame({}) # gait details 
    #TODO get gravity

#    date= '2017-10-12'    
#    dfsensor = vu.getVectorsBase(subject,date,
#                      rdf2,basedir,subdir,stime='00:00:00',etime='23:59:59')
#    dfsensor[['x','y','z']].plot()
    thrsd, thrrange = vs.onbodyThresholds[subject] # calculated in clinical         
    thrsd *= 1.00
    thrrange *= 1.00
    print('thresholds for subject,sd,range ',subject,thrsd,thrrange)

    
 
    smawindow = 50
    for i,fn in enumerate(adlfiles):
        # get file
       # fn = adlfiles.iloc[0].fn
        print('Subject:',subject,'processing file ',i,' of ',len(adlfiles),': ',fn)
        dfsensor = fh.getMatlabfile(fn,basedir,subdir,subject) 
        # dfsensor = returndf.copy()
        # get onbody periods 
        
        onbodywindow = 2 # minutes window size 
       #TODO standardise all times in seconds 
        onBodyPeriods, offBodyPeriods = vs.dissectData(dfsensor,onbodywindow,
                                               thrsd=thrsd,thrrange=thrrange,
                                               samplerate=samplerate)
        
        print('no. onbody periods for file ',fn,':',len(onBodyPeriods))
#        gravlist = vs.getGravity(dfsensor,offBodyPeriods,onbodywindow,samplerate)
#       this method produces results that seem too variable - using the value
#       calaculated inclinical tests instead 
#        gmean = np.nanmean(gravlist)
#        print('gravity',gmean)
#        gravity = gmean
        thisg = vs.gravity[subject]
        for periodnumber, period in enumerate(onBodyPeriods.values()):
            # get activity levels for each segment 
           # period = onBodyPeriods[2]
           # periodnumber = 1
            startperiod = int(period[0]*60*onbodywindow*samplerate)
            endperiod = int((period[1]+1)*60*onbodywindow*samplerate)
            print('on body period:',startperiod,endperiod,fn)       
            ondata = dfsensor.iloc[startperiod:endperiod,:].copy() 
            ar, arb, argr, ars, vacc, Vhacc, newdf = vs.getVectors(ondata,
                                                                   start=0,
                                                                   end=-1,
                                                                   gvalue=thisg)
            epoch = [newdf.datetime.values[0],newdf.datetime.values[-1]] #TODO these are np format - y?

            sma = vs.getSMA(arb)
            starttod = newdf.iloc[0].timestamp
            endtod = newdf.iloc[-1].timestamp
            print(datetime.utcfromtimestamp(starttod).strftime('%d/%m %H:%M:%S'))
            print(datetime.utcfromtimestamp(endtod).strftime('%d/%m %H:%M:%S'))

            t = np.linspace(starttod, endtod, len(sma))
            
            t2 = np.array([datetime.utcfromtimestamp(x) for x in t])

            actclasses=vs.getActivityClasses(sma,g=thisg) #... compute the activity classes 
            actdfs = pd.Series({'subject':subject,
                                'fn':fn,'sma':sma,'period': periodnumber,
                                'actclasses':actclasses,'t':t2})
            ADLactdf = ADLactdf.append(actdfs,ignore_index=True)
            datetime.utcfromtimestamp(starttod).strftime('%d/%m %H:%M:%S')
       
            if walksgait: # TODO refactor 
            # get walking bouts for each period of onbody time 
                statslen = 25
                walklog = vs.getWalks(vacc,actclasses,smawindow,epoch,samplerate=50,
                                      accelthov=1.4,statslen=statslen)
           
                for walk,data in walklog.items():
                    walkdfs = pd.Series({'subject':subject,
                                    'fn':fn,'period': periodnumber,
                                     'walk': walk,
                                    'startsec':data[0],
                                    'endsec':data[1],
                                    'nosteps':data[2],
                                    'actavg':data[3],
                                    't':data[4]
                                })
        
                    ADLwalkdf = ADLwalkdf.append(walkdfs,ignore_index=True)
              # get gait parameters for walks > 25 steps that have an sd < 0.02
                    for i,walksd in enumerate(data[5]):
                        
                        if i != 1: continue
                       
                        if walksd[2] < threshsd:
                           
                            # get gait features for this walk segment. the segment
                            # start/end is at walksd[0] & [1]
                            axyz = arb[walksd[0]:walksd[1]]
                            row = pd.Series({'walkid':str(walk)+'_' + str(i),
                                             'start':walksd[0], 'Method1_stepcount':statslen})           
    
                            rc, features =  vs.getGaitfeature(axyz,thisg,samplerate,row,plot_test=False)
     
                            # add gait features to adlgait df 
                            if rc != 1:
                                 print('bad rc from subwalk:',i)
                                # walkerrs.append(walk)
                            else:
                                print(features)
                                ADLgaitdf = pd.concat([ADLgaitdf,features])
                                print('walksd',walksd)
                                filepath = basedir + subdir + subject + '/ACT.pickle' 
    
    filepath = basedir + subdir + subject +'/ACT.pickle'
    pickle.dump(ADLactdf, open(filepath,'wb'))
    # save walks 
    if walksgait: 
        filepath = basedir + subdir + subject +'/walk.pickle'
        pickle.dump(ADLwalkdf, open(filepath,'wb'))
    
    return ADLactdf, ADLwalkdf,ADLgaitdf
            
       # we now hav an sma and act class for eaach on body segment for this id 
            
basedir = 'C:/Users/telferm/projects/vital/data/'
subdir = 'Batch1/'
subject = 'vah002'
global basedir
  
   # del rdf
    try:
        rdf
    except NameError:
        report, rdf = fh.createFilereport(basedir,subdir)

ADLactdf, ADLwalkdf,ADLgaitdf = manageActvitylevels(subject,rdf,basedir,subdir)
          
        #getInactiveSegments
       # save SMA  
      
        
    
#%% activity levels - group by time period 
def getsubjectSMAgrouped(subject,freq='H'):
# get sma levels from pickle  
    basedir = 'C:/Users/telferm/projects/vital/data/'
    subdir = 'Batch1/'
    #subject = 'vah006'
    subjectlevel ={}
    subjectlevel['vah006']= [1,2.5,5,7]
    subjectlevel['vah001']= [0.5,1.25,2.5,3.5]
    subjectlevel['vah002']= [1,2.5,5,7]
    subjectlevel['vah010']= [1,2.5,5,7]
    def calctime(x):
       plink =  x['stime'] + \
              timedelta(seconds=int(x['duration'])*60)
       return plink
    
    def addHour(x):
       plink =  x['stime'] + \
              timedelta(hours=1)
       return plink
    
    def classif(x,actlevels=[1,2.5,5,7]):
        if x < actlevels[0]: return 0 
        if (x > actlevels[0]) & (x <= actlevels[1]): return 1
        if (x > actlevels[1]) & (x <= actlevels[2]): return 2
        if (x > actlevels[2]) & (x <= actlevels[3]): return 3
        if (x > actlevels[3]): return 4
    filepath = basedir + subdir + subject + '/ACT.pickle' 
    actdf =  pickle.load(open(filepath,'rb'))
    # with dictionary of sma , for each plot the 5 minute smoothed intervals
    # actdf has one line per period of activity, so one day may be split across
    #  multiple files. 
    # amalgamate all sma files
    actdfall = pd.DataFrame({}) 
    
    for idx, ser in actdf.iterrows():
      #  if idx > 3: break
        combactdf = pd.DataFrame([ser.sma,ser.t,actdf.subject]).transpose()
        actdfall = pd.concat([actdfall,combactdf],axis=0)
    
    actdfall.info()
    
    actdfall.columns=['sma','t','subject']
    actdfall['numsma'] = pd.to_numeric(actdfall.sma)
    actdfall['numsma2'] = pd.to_numeric(actdfall.sma**0.5)
    # group into n minute intervals 
    ggg = actdfall.groupby(pd.Grouper(key='t',freq=freq)).mean()
    #ggg30min = actdfall.groupby(pd.Grouper(key='t',freq='30min')).mean()
    
    #fig, ax1 = plt.subplots()
    #ax1.plot(ggg.index, ggg.numsma)
    #ax1.fill_between(ggg.index,0, ggg.numsma)    
    # get the activities recorded 
    actdiarypath = basedir +'clinical/vitalhome_' + \
       subject + '_activity.csv'
    actentries = pd.read_csv(actdiarypath)
    #calculate start and enddates 
    
    actentries['stime'] = actentries['DATE'].map(lambda x: \
              datetime.fromtimestamp(x/1000))
    actentries['etime'] = actentries.apply(calctime,axis=1) 
    actentries['localstime'] = actentries.stime.apply(lambda x: x - timedelta(hours=1)) 
    actentries['localetime'] = actentries.etime.apply(lambda x: x - timedelta(hours=1))
    # map from text to numeric levels
    intensitymap = {'VERY_LIGHT':1,
                    'LIGHT':2,
                    'SOMEWHAT_HARD':3,
                    'HARD':4,
                    'VERY_HARD':5}
    actentries['intense'] = actentries.intensity.map(intensitymap)
    
    return ggg

smagroups = {}
subjects =['vah001','vah002','vah006','vah010']
weights =[91,110,69,79]

for subject in subjects:
    #subject = 'vah006'
    smag = getsubjectSMAgrouped(subject,'1H')
    smagroups[subject] = smag.copy()
  
sdatev006 = smagroups['vah006'].index.min() 
fig1, ax1 = plt.subplots()
for subject,dfs in smagroups.items():
#    print(subject)
#    subjectix = subjects.index(subject)
#    wt = weights[subjectix]
#    #move start date to vah006 (earliest)
#    # calc vah start date
#    
#   
#    #sdate = type(startdate)
#    # get dif in days for this subject 
#    sdate =  smagroups[subject].index.min()
#    diff = sdate - sdatev006
#    diff = pd.Timedelta(days=diff.days)
#    print('subject',subject,'days',diff)
#    # substract from all dates to create new date field 
#    dfs['rebasedstart']=dfs.index - diff 
#    
#    dfs['smaroot']=dfs.numsma**0.5 # de power the sma to energy 
#    dfs['energy']=dfs.numsma*wt # de power the sma to energy 
    ax1.plot(dfs.rebasedstart,dfs.numsma,label=subject)
    fig1.suptitle('SMA values - average over 1 Day')
plt.show()
ax1.legend()
# end time = start time + 60 * duration in seconds 

# set level according to classif    
#ggg['level'] = ggg.numsma.map(lambda x: classif(x,subjectlevel[subject]))
#minsma = ggg.numsma.min() # to translate to 0
#maxsma = ggg.numsma.max() - ggg.numsma.min() # will be max after translation
#
#ggg['scaledsma'] = ggg.numsma.map(lambda x: (x-minsma)*(5/maxsma))    
# 
#ggg['scaledsma'].max()
 


#%%plot sma against activity levels 
def make_patch_spines_invisible(ax):
    '''trying to get the sma real values on the x axis and failing '''
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return
        
        
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

for idx,x in actentries.iterrows():
    ax2.plot_date([x['localstime'],x['localetime']],[x['intense'],x['intense']]
    ,'ro-')
ax2.plot_date([],[] ,'ro-',label='Diary Entry')

ax1.plot(ggg.index, ggg.level,label='SMA derived levels')
plotreal = False

if plotreal:
    ax3 = ax1.twinx()    
    realsmaline = ax3.plot(ggg30min.index,ggg30min.numsma,'r-',
                           alpha=0.5,label = 'SMA real value 2')
    #pos1 = ax3.get_position() # get the original position 
    #pos2 = [pos1.x0+0.1, pos1.y0,  pos1.width, pos1.height]
    ax3.spines["left"].set_position(("axes", -0.2)) # green one
    make_patch_spines_invisible(ax1) 
    ax3.spines["left"].set_visible(True)
    
    ax3.yaxis.set_label_position('left')
    ax3.yaxis.set_ticks_position('left')
    ax3.set_ylabel('real sma values')
    ax3.set_ylim(0,int(round(maxsma,0)))
    ax3.set_yticks([0,1,2,3,4,5,6])
    ax3.legend(loc=9)

#r = realsmaline.pop()
#r.remove()

#ax1.fill_between(ggg.index,0, ggg.numsma2) 
ax1.set_ylabel('SMA derived levels')
ax1.set_xlabel('Date Time')
ax2.set_ylabel('Self ascribed activity labels')

ax2.set_ylim(0,5)

#ax1.set_ylim(0,30)
ax2.set_yticks([0,1,2,3,4,5])
ax1.set_yticks([0,1,2,3,4,5])
#ax1.set_zorder(1)
#ax2.set_zorder(2)
#ax3.set_zorder(3)
ax1.patch.set_alpha(0.1)
ax2.patch.set_alpha(0.1)
ax1.legend(loc=2)
ax2.legend()
hourFmt = mdates.DateFormatter('%H:%M')
daysFmt = mdates.DateFormatter('%d/%m')
hours = mdates.HourLocator((8,12,18) )  # every hour
days = mdates.DayLocator(interval=1)   # every day
ax1.xaxis.set_major_formatter(daysFmt)
ax1.xaxis.set_major_locator(days) 
ax1.xaxis.set_minor_formatter(hourFmt)
ax1.xaxis.set_minor_locator(hours)  
fig.autofmt_xdate()
ax1.xaxis.set_tick_params(which='minor',labelrotation=30,labelsize=6)
for label in ax1.get_xminorticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment("right")
ax2.set_yticklabels(['']+[x for x in intensitymap.keys()],rotation=60,va='bottom')
fig.suptitle(subject + ': using levels '+ str(subjectlevel[subject]))
plt.show()
#%% calculate no. minutes in each bucket : 

    
ggg['level'] = ggg.numsma.map(classif)
ax1.plot(ggg.index, ggg.level)
ax1.set_ylim(0,5)
ax1.set_ylabel('SMA derived levels')
ax1.legend()

# plot overall freq

tt = ggg.level.value_counts().values

type(tt) 
tt = tt[1:] 
tt = tt*5    
minutesdiary = (20,280,140,0)
ind = np.arange(4)  # the x locations for the groups
width = 0.35 
fig = plt.figure()
ax = fig.add_subplot(111)
fig,(ax,ax2) = plt.subplots(2, 1, sharex=True)

rects1 = ax.bar(ind, tt, width, color='r')
rects1 = ax2.bar(ind, tt, width, color='r')

rects2 = ax.bar(ind+width, minutesdiary, width, color='y')
rects2 = ax2.bar(ind+width, minutesdiary, width, color='y')

# add some
ax.set_ylabel('Minutes')
ax2.set_xlabel('Activity Level')
ax.set_xticks([])
#ax2.xaxis.set_major_formatter(plt.NullFormatter())
ax.set_title('Minutes Activity by Sensor vs Diary')
#ax2.set_xticks(ind+(width/2))
ax2.set_xticks([0,1,2,3,4])
ax2.set_xticklabels( ('1', '2', '3', '4') )
ax2.set_ylim(0,500)
ax.set_ylim(1500,2000)
# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.yaxis.tick_left()
ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.yaxis.tick_left()

ax.legend( (rects1[0], rects2[0]), ('SMA', 'Diary') )

plt.show()
