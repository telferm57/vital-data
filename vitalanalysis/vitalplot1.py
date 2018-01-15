# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 05:40:15 2017

@author: telferm
"""


#%% plotting helpers
import matplotlib.pyplot as plt
from matplotlib import patches 

import pandas as pd
from datetime import timedelta
from datetime import datetime
import matplotlib.dates as mdates
import numpy as np

def addGrid():
    from matplotlib.ticker import AutoMinorLocator
    ax1 = plt.gca()
    ax1.minorticks_on()
    minor_locator = AutoMinorLocator(10)
    ax1.xaxis.set_minor_locator(minor_locator)
    plt.grid(which='both')
   # ax.legend()

#ax.set_xticks(tst.values)
def pltAnnotations(ax,annotations):
    annottimes=mdates.date2num(annotations.synched.tolist())
    ax2 = ax.twiny()
    #ax.set_xticks(annottimes)
    #ax.set_xticklabels(annottimes,rotation = 45)
    ax2.set_xticks(annottimes)
    ax2.set_xticklabels(annotations.notes2.values,rotation=80,ha='left')
    #ax.set_xlim(annottimes[0],annottimes[-1])
    #ax2.set_xlim(tst[0],tst[-1])
    ax2.set_xlim(annottimes[0],annottimes[-1])
    
   # timeFmt = mdates.DateFormatter('%M:%S')
   # ax.xaxis.set_major_formatter(timeFmt)
    out =  ax2.plot() #,ax.plot()
    return ax2,annottimes
   
def plotVector(v,epoch,window):
    ''' create x axis for vector v , starting at start of epoch 
    assuming each element is window seconds from last. place each
    value in the middle of the window. window defined in seconds '''
    # so for SMA we have window of 1 second 
    window = 1
    td = timedelta(seconds = window)
    base = datetime(2000, 1, 1)
    timearray = np.arange(epoch[0]+np.timedelta64(seconds=window/2), 
                          epoch[1],np.timedelta64(window,'s'), 
                          dtype='datetime64')
    # convert to datetime 
    
    return pd.to_datetime(timearray)

def plotHeatAct():
    ''' heatmap for activity level - this has never been run as a function and needs tidying up
    considerably ! '''
    actdata =  pd.read_csv('activiyADL006.csv',parse_dates=True)

    acttt = actdata[['date','time','sma']].copy()


   # remove invalid dates from date column 
    def nainvaliddates(x):
        try:
           cnv_dt = datetime.strptime(str(x), '%Y-%m-%d')
        except ValueError:
           cnv_dt = np.nan
        return cnv_dt


    acttt['date2'] = actdata.date.map(nainvaliddates)
   
    acttt.shape

    acttt.dropna(subset=['date2'],inplace = True)

   

    import seaborn as sns
   # acttt.date.iloc[1:2]
    from datetime import time

  
    sdate = datetime(2017,9,21,hour=4)
    sdate = time(hour=5,minute=0,second=0)
  
    acttt['dtime'] = acttt.time.apply(lambda x : datetime.strptime(x,'%H:%M:%S').time())
    acttt['ddatatime'] = acttt.apply(lambda x: datetime.strptime(x['date'] +          \
                               'T' + x['time'],    \
                               '%Y-%m-%dT%H:%M:%S'),axis=1)
    acttt.info()
   
    actttshort = acttt[acttt.dtime > sdate].copy()
    actttshort.info()
   # need to group the reading in 1 minute intervals 
    ggg = actttshort.groupby(pd.Grouper(key='ddatatime',freq='1Min')).mean()
    ggg.info()
    ggg['date'] = ggg.index.map(lambda x: datetime.date(x).strftime('%Y-%m-%d'))
    ggg['time'] = ggg.index.map(lambda x: datetime.time(x).strftime('%H:%M:%S'))

    actpivot = pd.pivot_table(ggg,values='sma',columns='date', index = ggg.time, 
    dropna=True)

    sns.heatmap(actpivot,yticklabels=10)
    
    return

def plotWalkHeat(walkdf):
    
    tt = walkdf.copy()
    tt = tt[tt.nosteps>3].copy()
    tt['duration'] = tt.endsec - tt.startsec

    ggg = tt.groupby(pd.Grouper(key='t',freq='H')).sum()
    #create date and time column
    ggg['date'] = ggg.index.map(lambda x: datetime.date(x).strftime('%Y-%m-%d'))
    ggg['time'] = ggg.index.map(lambda x: datetime.time(x).strftime('%H:%M:%S'))
    walkpivot = pd.pivot_table(ggg,values='duration',columns='date', index = ggg.time, 
    dropna=True)

    ax = sns.heatmap(walkpivot,yticklabels=6,linewidths=0.5,)
    ax.set_title('VAH010: seconds of walks per hour per day')
    
    ax.set_yticklabels(ax.get_yticklabels(),rotation=0,fontsize=8)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
    ax.invert_yaxis()
    return                    


def plotWalkTimes(ax, walklog,epoch,rate=50, patchy=1,param_dict={}):
    from datetime import timedelta
    import matplotlib.dates as mdates 
    ''' given a walklog object, that records all walks over a predfined sensor
    epoch, and a set of axes, draw the blocks of time on the axes 
    corresponding to walks > 2 steps. walks < 2 steps may be something else 
    
    input: walklog - dict of walks that have start& end  time (seconds) and
    number of steps per walk 
        epoch : list of datetime object for the start and end of the epoch under examination
        
    #TODO may be better to have these indexed by time rather than sample numbers 
    which probably means changing the walklog to have start time, endtime rather 
    than seconds since the beginning 
    '''
   # walklog=walks.copy()
    # get min max for walk speeds to get rgb limits 
    ll = [(v[1]-v[0])/v[2]  for k,v in walklog.items() if v[2] > 2]
    # scale min-max for g=rgb 0.25 - 0.75
    vrange = max(ll) - min(ll)
    vmin = min(ll)
    ll2 = ((ll-vmin)/vrange)*0.5+0.25
    
   
    
    stime = pd.to_datetime(epoch[0])
    
    for walk,detail in walklog.items():
        if detail[2]> 2:  # walklength(no. steps)  > 2 
            print('walklog ',walk)
            stt = stime + timedelta(seconds=detail[0])
            ent = stime + timedelta(seconds=detail[1])
            x = mdates.date2num(stt)
            x2 = mdates.date2num(ent)
           
            w= (detail[1]-detail[0]) # length walk in sec
            h=patchy*2
            #x = 20000
            y = -patchy
           # w=20000
            steptime= w/(detail[2]) # in seconds
            #print(steptime,type(steptime))
            wt = x2 - x # width in matplotlib time format  
            rgbrval = ((steptime-vmin)/vrange)*0.8+0.2
          #  print(rgbrval,steptime,type(steptime))

            ax.add_patch(patches.Rectangle((x,y),wt,h
                                           ,fc=(rgbrval,0,0)
                                           ))
            

            ax.text(x,patchy,str(detail[2])+'('+str(walk)+')',
                    color='lightgreen',va='top')
            ax.text(x,-patchy,'{0:1.2f}'.format(steptime),color='lightblue',
                    rotation=90,va='bottom')        
    out = ax.plot()
    return out

def plotPTTimes(ax,PTlog,PTdetail,epoch,rate=50,offset=0,param_dict={}):
    '''  plot pt times on axes 
    note that the x axis passed in has tod. the PTlog is in samples, so need to convert samples to timedelta 
    in seconds and add to epoch start '''
   
    h=20
    y=-10

    stime = pd.to_datetime(epoch[0])
    # PTs last between 0.5 and 4 seconds .. so throw away anything 
    # 1/2 the rate and 4 * rate- although are PTs contigous ? they are in FTSTS 
    #cut down dict 
    #for PT in (y for y in PTlog.values() if (y[2] >= 20)):
    for PT in [[k]+y for k,y in PTlog.items() if (y[2] >= 20)]:
       # PT = PTlog[272] # 272, 341 
       #print(PT)
        x,x2 = PT[1]/rate,PT[2]/rate
        stt = stime + timedelta(seconds=x)
        ent = stime + timedelta(seconds=x2)
        x = mdates.date2num(stt)
        x2 = mdates.date2num(ent)
    
        # duration of transition
        w= x2-x
        # get transition type from ptdetail 
          
        PTtype = PTdetail.get(PT[0],0) 
        if PTtype:
            PTtype = PTtype[0]
            fc = 'b'
        else:
            PTtype = -1 # indicates invalid according  to PTdetail 
            fc = 'g'
        ax.add_patch(patches.Rectangle((x,y),w,h,
                                           fc=fc,alpha=0.2))
        #add classification 1 or 2 
      
        
    
        
    
       # PTtype = PTdetail.setdefault(1,-1) 

        ax.text(x,y+h,str(PTtype),color='black',
                    rotation=0,va='bottom')   
        
        
    out = ax.plot
    return out
#plt.show()

def plotPTdetail(ax,PTdetail,epoch,rate=50,param_dict={}):
    ''' plot PT details on axes - intended to overlay the Vacc and sintheta 
    plots 
    
    PT detail consists of the classification of the PT, 3 x,y 
    co-ordinates of the sintheta peak and the vert velocity trough and peak 
    and the ratio of peak to trough. The ratio is used to classify the 
    transition '''
    # get first set of points
    # PT detail is in samples, not time - need to convert to
    # display on time axis 
    PTclass = [x[0] for x in PTdetail.values()]
    sinth = [x[1] for x in PTdetail.values()]
    PTratio = [x[4] for x in PTdetail.values()]

    vvpeaks = [x[2] for x in PTdetail.values()]
    vvvalleys = [x[3] for x in PTdetail.values()]
    sinthpts = [[x[0] for x in sinth],[x[1] for x in sinth]]
    vvpeakspts = [[x[0] for x in vvpeaks],[x[1] for x in vvpeaks]]
    vvvalleyspts = [[x[0] for x in vvvalleys],[x[1] for x in vvvalleys]]
    ax.plot(sinthpts[0],sinthpts[1],'o', mfc=None, mew=2, ms=8,
             label='Sintheta peaks')
    ax.plot(vvpeakspts[0],vvpeakspts[1],'+', mfc=None, mew=2, ms=8,
             label='Vert Velocity peaks')
    ax.plot(vvvalleyspts[0],vvvalleyspts[1],'+', mfc=None, mew=2, ms=8,
             label='Vert Velocity valleys')

    
    for i,p in enumerate(PTclass):
        tb = ax.text(sinth[i][0],sinth[i][1]*1.1,str(p)+
                     ' '+'{0:1.2f}'.format(PTratio[i]))
        tb.set_bbox(dict(facecolor='red', alpha=0.5, edgecolor='red'))
                

    
    out = ax.plot()

    return out

def plotInact(ax, iasegs,rate=50,param_dict={}):
    '''  plot pt times on axes  '''
    h=4
    
    for ia in iasegs:
        x,x2 = (ia[0]-1)*rate,(ia[1])*rate
            
        w= x2-x + 1
        print(x,w,x)
        y=-2
        ax.add_patch(patches.Rectangle((x,y),w,h,
                                           fc='y',alpha=0.3))
    out = ax.plot
    return out

def plotPressure(ax, pressure,rate=25,param_dict={}):
    plt.plot(periodpressure.values,label='pressure')


