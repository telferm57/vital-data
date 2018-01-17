# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:12:47 2017

@author: telferm
"""
#from vitalanalysis import filehandler as fh
#from vitalanalysis import vitalsensor as vs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.signal import medfilt
from scipy.integrate import cumtrapz 
import imp 
#imp.reload(mhealthx.extractors.pyGait)
from math import floor#import vitalanalysis.cfg as cfg
from vitalanalysis import fh
#import vitalanalysis.filehandler as fh
import pandas as pd
from datetime import timedelta
from mhealthx.extract import run_pyGait
import mhealthx.extractors.pyGait
import mhealthx.signals
from mhealthx.extractors.pyGait import project_walk_direction_preheel, walk_direction_preheel
global onbodyThresholds
onbodyThresholds = {} 
global gravity
gravity = {} 
#%% signal processing stuff
from scipy.signal import butter, lfilter, welch, periodogram, find_peaks_cwt, filtfilt
def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_bandpass_filterf(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band',analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=4):

    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=5):

    nyq = 0.5 * fs
    cutf = cutoff / nyq
    b, a = butter(order, cutf, btype='low',analog=False)
    return b, a

def lpFilter(x,y,z,cutoff,samplerate,order=4):

    xfilt = butter_lowpass_filter(x,cutoff,samplerate,order)
    yfilt = butter_lowpass_filter(y,cutoff,samplerate,order)
    zfilt = butter_lowpass_filter(z,cutoff,samplerate,order)
    return xfilt,yfilt,zfilt

def bpFilter(x,y,z,lowcut,highcut,samplerate,order=4):
    xfilt = butter_bandpass_filter(x,lowcut,highcut,samplerate,order)
    yfilt = butter_bandpass_filter(y,lowcut,highcut,samplerate,order)
    zfilt = butter_bandpass_filter(z,lowcut,highcut,samplerate,order)
    return xfilt,yfilt,zfilt


#def magnitude(ar): # assuming x, y, z array 
#''' given an array with 3 coumns or rows, return 1d array of magnitude 
#along the longest axis '''
#      
#    ar2 = np.square(ar)
#    m2 = ar2[0] + ar2[1] + ar2[2]
#    m = np.sqrt(m2)
#    return m

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

#%% Identify segments of interest 
def getTestData(fn,nrows=1e4,skiprows=0):
    import pandas as pd
    path='C:/Users/telferm/python/projects/walkingsignals/walkingsignals/data/smm/pilot/' 
    smm = pd.read_csv(path + fn,header=None,skiprows=skiprows,nrows=nrows)
    return smm 


import numpy as np
def getRange(ar):
    ''' calculate sd and range for 3 axes in a resting state.
    typically to find the offbody segments, I use a medfilt of 3 '''
    return ar.max()-ar.min()

def getMonitorData(subject,rdf,basedir,adl=0):
    ''' only retrieves the one file for the visit ! not ADL 
    TODO incorporate into other routines  - is this a filehandler routine ? 
    '''
    print('retreiving clinical monitor data for ',subject)
    fn = rdf[(rdf.subject==subject) & (rdf.adl == adl)].fn.values
    #should only be one 
    fn = fn[0]
    subdir = 'Batch1/'
    #path = basedir + subdir
    
    #get matlabdata 
    dfsensor = fh.getMatlabfile(fn,basedir,subdir,subject)
    
    return dfsensor

def getMonitorThresholds(subject,rdf,basedir,startsample,endsample,dfsensor=False):  
    '''
    Calculate the range and std for the monitor at rest. This is then used in 
    other functions to calculate the on body/off body time segments 
    Calculates gravity using getGravity function
    
    input: 
    subject, rdf (sensor file names, subjects, dates) 
    startsample, endsample - sample number of start/end of inactive period (i.e.
    when the sensor is off the body )
    dfsensor - dataframe containing sensor output in which there is an inactve
                period. If not a dataframe, will the clinical data using
                getMonitordata 
    output:
    sets gobals gravity, onbodyThresholds
        '''
#TODO change to get file report and create if necessary
    global onbodyThresholds
    global gravity 
    if not isinstance(dfsensor,pd.DataFrame): 
        print('acquiring monitor thresholds for subject:',subject)
        dfsensor = getMonitorData(subject,rdf,basedir)
    #dfsensor = senscomb.copy()
    # find low segments to assess thresholds for off body 
    #dfsensor[['x','y','z']].plot()
    samplerateacc  =len(dfsensor)/(dfsensor.timestamp.max()-dfsensor.timestamp.min())
    #samplerate = np.around(samplerate,decimals=1)
    samplerate = int(round(samplerateacc,0))
    # get off body parameters - max sd and max range while off body 
    # oding this by eye as we have no defined time when the sensor is off the body
    #
    restdata = dfsensor[['x','y','z']][startsample:endsample]
    
    thrsd, thrrange = getNoise(restdata,subject) # get thresholds for testing sensor off body
    # get gravity value fr this segment 
    onbodyThresholds[subject] = [thrsd, thrrange]
    
    gravity[subject]=getGravity(restdata,False,0,samplerate)[0] 
    print('gravity for ',subject,' calculated as ',gravity[subject])
    return  gravity[subject], thrsd, thrrange, samplerateacc
    
def testOffBody(data, thrsd=0.1,thrrange = 0.5):
   
    ''' given input accelerometer data, return whether device was off body or not
     #      https://doi.org/10.1371/journal.pone.0022922 Van Hees 2011
     
     2 axes have a std deviation < 0.5; OR 2 axes have a range < 0.8 
     
    " A block was classified as non-wear time if the standard deviation was
     less than 3.0 mg (1 mg = 0.00981 m·s−2) for at least two out of the three
     axes or if the value range, for at least two out of three axes, was less 
     than 50 mg"
     
     3mg = 0.03 m/s 50 mg = 0.5 m/s 
     
    n.b. for smm, the sd and range can exceed these values when the device is left on 
    a table, for example.
    
    e.g. for vah006, sd averages 0.094, range averages 0.58 (both values calculate 
    after medfilt3)
    
        
     returns 0: on body
             1: off body 
    Have now changed this to use getNoise to calculate reasonable 
    thresholds for each device (the data for the getnoise must be identified manually)
    '''
    def getRange(ar):
        return ar.max()-ar.min()
    
    res = 0 # it IS on body unless ...
    axsdtot = 0; axrangetot = 0 
    #print(data.columns)
    sds = {}
    ranges = {}
    for col in data.columns:
       
        sdaxis = np.std(data[col].values)
        sds[col] = sdaxis
        if sdaxis < thrsd: 
            #print(col,' sd gt 0.03 :',sdaxis)    
            axsdtot += 1
        #print(col, sdaxis, getRange(data[col].values))           
        if getRange(data[col].values) < thrrange:  
            #print(col,' sd gt 0.03 :',sdaxis)    
            axrangetot += 1
    if (axsdtot > 1) |  (axrangetot > 1) : res = 1         
    #print(col,' no.  axes exceeded 1, sd, range :',axsdtot,axrangetot)    
    return res,sds

def getOffBody(ws,data,rate,thrsd=0.1,thrrange=0.5):
    ''' for a given window size, return a dict of windows indicating 
    whether device on (0) or off (1)  body '''
#    data=tt.copy()
#    ws = 2 # minutes
#    rate=50
    offbody = {}
    allsds = {}
    samplesize = int(round(ws*rate*60,0)) # 6000 per 2  minutes 
    print('no. 2 minute samples:',floor(len(data)/samplesize))
    for sampleno in np.arange(floor(len(data)/samplesize)):
        sstart = sampleno*samplesize
        send = sstart+samplesize
        sample = data[sstart:send].copy()
        offbody[sampleno],sds = testOffBody(sample,thrsd,thrrange)
        for col in sds.keys():
            allsds.setdefault(col,[]).append(sds[col])
        
            
    return offbody,allsds  

def dissectData(data,onbodywindow=2,thrsd=0.1,thrrange=0.5,samplerate=50):
    ''' for a given datastream, identify on-body segments 
    [and maybe pickle them
    input: pandas dataframe of data with accelerometer in columns named x,y,z
            window size in minutes 
    
    returns number of periods dict '''
    
    tt= data.copy()
    #   tt= returndf.copy()[30000:-60000]
    tt = tt[['x','y','z']]
   # tt.columns = ['x','y','z'] 
    if tt.isnull().sum().sum() != 0:
        print("nulls in data - I can't handle that ")
    ar = tt[['x','y','z']].values 
    smoothedx = medfilt(ar[:,0],3)   
    smoothedy = medfilt(ar[:,1],3)   
    smoothedz = medfilt(ar[:,2],3) 
    
       
    #tt['sx','sy','sz'] =smoothed[:,0:2] 
    tt['sx'] =smoothedx
    tt['sy'] =smoothedy 
    tt['sz'] =smoothedz         
   # tt.plot()    
    #np.std(tt.z.values)
    #tt.describe()
    tt.drop(['x','y','z'],axis=1,inplace=True)
    #tt.plot()
    #samplerate=50

    offBodySegments,allsds = getOffBody(onbodywindow,tt,samplerate,
                                        thrsd,thrrange) # x minute mindow, 50 samples/second 
   #plt.plot(np.asarray(list(offBodySegments.values())))
    #plt.show()
    #offBodySegments = offBodyPeriods.copy()
    offBodyar = np.asarray(list(offBodySegments.values()))
    # get periods of segments from ' off body' dict 
    
    #only interested in > 10 minutes of wear time/offtime, so 5 consecutive zeros or
    # TODO 5 consecutive ones 
    #offBodySegments = offBodyPeriods.copy()
   
    runsize = 5
    periodcounter = 0 
    periods = {} # generate dict of periods of on body - runs of 0's 
    switch = np.sign(offBodyar[:runsize].sum()) # initial state 0 if 0 , 1 if > 0 
    if switch == 0: startrun=0
    for i in np.arange(len(offBodyar)-runsize):
        runsum = offBodyar[i:i+runsize].sum()
        # at beginning , record beginning
        oldswitch = switch 
        switch =  np.sign(runsum) # either 0 or > 0 
        if oldswitch != switch: # state transition
            if switch == 0: # start of new period
                print('start new period ',i)
                startrun = i # all elements must be zero 
            else: # end of period switch from 0 to 1 - on to off body
                endrun = i+ runsize-2 # the first runsize-1 values are actualy 0 
                print('end period ',endrun)
                periods[periodcounter] = [startrun,endrun]
                periodcounter += 1
        if (i+runsize == len(offBodyar)-1) & (switch==0): # end of array while in 0
            print('end of array processing',i+runsize,len(offBodyar)-1,switch) 
            endrun=i+runsize
            periods[periodcounter] = [startrun,endrun]
            
            
            

    ws=2 # minutes
    samplesize = int(ws*samplerate*60)
  
    def picklePeriods(data, periods, samplesize,fn):
        import pickle
        import os
        epoch={}
        #os.getcwd()
        for period in periods:
            length = periods[period][1]-periods[period][0] + 1
            length *= samplesize
            start = samplesize*periods[period][0]
            epoch[period] = data[start:start+length] 
            pickle.dump(epoch[period], open(fn + '_' + str(period),'wb'))
            
        return
                
    #fn='userid-date-device'            
    #picklePeriods(data, periods, samplesize,fn)
   # import pdb
    #pdb.set_trace()
    
    return periods, offBodySegments

def getNoise(restdata,subject):
    ''' return the mean std and range for 2  minute windows for monitor at rest.
    values are median filtered prior, as they are in the rest of this app''' 
   
    ar = restdata[['x','y','z']].values 
    ar.shape
    samplerate=50
    sx = medfilt(ar[:,0],3)   
    sy = medfilt(ar[:,1],3)   
    sz = medfilt(ar[:,2],3) 
    sxyz=np.asarray([sx,sy,sz]).T
    
    sd = []
    
    rangevalues = []
    for (start, end) in windows(sxyz,2*60*samplerate):
        win = sxyz[start:end,:]
        windf = pd.DataFrame(win)
        sd.append(list(windf.std().values))
        rangevalues.append(list((windf.max()-windf.min()).values))
    windf = pd.DataFrame(sd)
    thrsd = windf.mean().mean()
    windf = pd.DataFrame(rangevalues)
    thrrange = windf.mean().mean()
    return thrsd, thrrange

def getGravity(data,offBodyPeriods,onbodywindow,samplerate):
    ''' gravity reading will vary for each device.. this routine will calculate
    what the particular device is measuing as 1g. that value will be used 
    throughout to calculate thresholds.
    
    Eitheroff body periods can be used for this, or a specific portion of the
    signal can be passed in 
    
    '''
    gravlist = []
    ddd = data[['x','y','z']].values
   
    if offBodyPeriods:
        for key, val in offBodyPeriods.items():
            if val == 1:
                start = int(key*samplerate*onbodywindow*60)
                end = int(start + samplerate*onbodywindow*60 -1)
                ar = ddd[start:end]
                gmag = magnitude(ar.transpose())
                gest = np.mean(gmag)
                gravlist.append(gest)
    else:
        gmag = magnitude(ddd.transpose())
        gest = np.mean(gmag)
        gravlist.append(gest)
                          
    
    return gravlist

def getVectors(data,start=0,end=-1,gvalue=9.81):
    ''' from accelerometer data passed in as dataframe  , extract
    1. gravity vector 
    2. body vector
    3. vertical acceleration 
    4. realign original vectors with lag induced in the claculated vectors 
        by the butterworth filters 
    5. realigned horizontal vectors 
    
    input: dataframe of x, y , z ; start pos; end pos 
  
    '''
    #data=walk15.copy()
    from scipy.signal import medfilt

    #gvalue=9.82  # note g can change over time dependent on accelerometer 
    #gvalue=9.74 
    tt = data[start:end].copy()
    tt = data[['x','y','z']] 
    tt.isnull().sum()
    ar = tt[['x','y','z']].values 
    sx = medfilt(ar[:,0],3)   
    sy = medfilt(ar[:,1],3)   
    sz = medfilt(ar[:,2],3)   
    # now butterworth with recommended filter params : 
    # filter order 2, cut off frequency 1.6Hz delay 0.141 seconds 
    xg,yg,zg = lpFilter(sx,sy,sz,1.6,50,2)
    # delay signal by 0.141 s = 7 slots - chop 7 off beginning of true, 7 off the end of filtered
    xg = xg[:-7]
    yg = yg[:-7]
    zg = zg[:-7]
    sx = sx[7:]
    sy = sy[7:]
    sz = sz[7:]
    # do the same with the original data 
    ar = ar[7:,:]
    #and adjust the dataframe 
    newdf = data.iloc[7:,:].copy()
    # calculate inertial signals - i.e. body movements, by subtracting gravity 
    xb = sx - xg
    yb = sy - yg
    zb = sz - zg
    arb = np.asarray([xb,yb,zb]).transpose() 
    argr = np.asarray([xg,yg,zg]).transpose() 
    ars = np.asarray([sx,sy,sz]).transpose()
    # normalise and get get dot product  for vertical acc 
    normg = (xg**2+yg**2+zg**2)**0.5 # equiv to gdotg
    vacc = np.array([np.dot(a,b) for a,b in zip(argr,ars)])
    vacc = vacc/normg 
    vacc = vacc-gvalue   # just to be explicit - body only 
    # calculate in a different way using body acc only 
    # this way introduces a variation in G that seems to destabilise the
    # signal slightly - it seems more valid to use the fixed gravity constant
    # calculated for the device 
#    vacc2 = np.array([np.dot(a,b) for a,b in zip(arb,argr)]) 
#    vacc2 = vacc2/normg
#    vacc2 = np.array([vacc2]).transpose()
#    #vacc2.transpose().shape
#   # argr.shape
    vacc2d = np.array([vacc]).transpose()
    Vvacc = vacc2d*argr # vertical acceleration expressed as a vector 
    Vhacc = arb - Vvacc

    return ar, arb, argr, ars, vacc, Vhacc, newdf



def windows(ar, size=50,overlap=0):
    start = 0
    while start <len(ar):
        yield start, start + size
        start += size-overlap


def getSMA(data,winsize=50,gvalue=9.81):
    #TODO convert winsize to seconds, use computed samplerate 
    ''' calculate the signal magnitude area 
    
   #input - standard sensor dataframe - extract  input 3d array shape N,3 
    input: body acceleration values (i.e. no gravity included )
    output : array of SMA values for each window of windowsize samples 
    '''
    # axyz = data[['x','y','z']].values what ? 
    from scipy.integrate import cumtrapz 
    integrals = np.array([])
    integralsc = np.array([])
    
    for (start, end) in windows(data,winsize):
        
        # if start%1000 == 0: print(start,end) # windows look good 
        # start,end =2000,4000
        win = data[start:end,:]
          # integral = cumtrapz(win,axis=0) - also checking if cumsum makes a difference
        integral = cumtrapz(abs(win),axis=0)       
        integrals = np.append(integrals,sum(integral[-1])/(winsize))
        
        integralc = np.cumsum(abs(win),axis=0) # todo eval simpsons 
        integralsc = np.append(integralsc,sum(integralc[-1])/(winsize))
        
    return integrals # from plot, cumtrapz and cumsum almost identical 

def getInactiveSegments(sma,minlength,thresh=1.5):
    
    ''' return list of inative segments .. timescale is unknown
     to this routine 
     
     output : list of inactive segments greater than minlength long.
     segment counting starts at 1 '''
     #inasegs = {}
     # get indexes of segs greater than thresh 
    inact = np.less(sma,thresh)
    inactstart = 0 
    inactsegs = []
     
    for i,val in enumerate(inact,start=1):
        if val: # less than active thresh
            if inactstart:
                continue
            else: 
                inactstart = i 
        else: # active - > active thresh 
            if inactstart:
                if i-inactstart > minlength-1:
                    inactsegs.append([inactstart,i-1])
                inactstart = 0 
    
    return inactsegs

def getActivityClasses(sma,threshlo=0.16,threshhi=0.7,threshfall=2.1,g=9.81):
    
    ''' return list of activity level of segments .. timescale is unknown
     to this routine 
     see Davide Curone et al. 2010 for thresholds and calculations 
     output : list of activity levels of  segments.
     0- low level such as sitting down 
     1- medium - walking , light gardening
     2 - high - walking intensely, running etc 
     3 - very high - abnrmaly activity such as fall, taking device off and throwing 
     segment counting starts at 1 '''
     #inasegs = {}
     # get indexes of segs greater than thresh 
     
    fall = np.greater(sma,threshfall*g)
    hi = np.greater(sma,threshhi*g)
    lo = np.greater(sma,threshlo*g)
    #actstart = 0 
    #actclasses = []
   # reset segments < minlength to 0 , for each threshold  - smoothing  
# trouble with this approach is if there are lots of < minlength alternations , they will be set to zero
# so going for smoothing instead ..  of the input    
#    for thresh in [lo,fall,hi]:
#       
#       for i,val in enumerate(thresh,start=1):
#           if val: # > thresh
#               if actstart:
#                   continue
#               else: 
#                   actstart = i 
#           else: # active - > active thresh 
#               if actstart:
#                   if i-actstart+1 < minlength: # if not minlength, zeroise from actstart to i 
#                       thresh[actstart-1:i] = [0]*(i-actstart)
#                       actstart = 0 
 
    actclasses = lo*1+hi*1+fall*1
    # actclasses = 
    
    return actclasses
    
    
 # test - minlength 2    

#%% Walking related functions 

def initValues():
    global potentialBreach,breachState,breachStart,walkState,walkCount, \
     walklog,  \
     arb,       \
     stepCount,     \
     walkCount,     \
     blipSize,      \
     stepGapMin,    \
     walkGapMin,    \
     accelTh,       \
     accelThhi,      \
     maxBreachLen,  \
     minBreachLen,  \
     subjects,  \
     basedir,  \
     subdir,  \
     step
    # all times in seconds - converted by sectopoints
    global walkGapMax
    global gravity
    gravity = {}
    global onbodyThresholds
    onbodyThresholds = {}
    global stepGapMins, walkGapMaxs, maxBreachLens, minBreachLens
    subjects =['vah001','vah002','vah006','vah010']
    basedir = 'C:/Users/telferm/projects/vital/data/'
    subdir = 'Batch1/'
    potentialBreach = 0 # don't think we need this
    breachState = 0 
    breachStart=0
    walkState = 0 
    walkCount = 0
    walklog = {}
    stepCount = 0 # within walk
    walkCount = 0 
    blipSize = 2 # number of readings that constitute a blip 
    stepGapMin = 4 # min number of readings between breaches to delineate step 
    stepGapMins = 0.4 # min number of readings between breaches to delineate step 
    walkGapMax = 25 # max  gap between steps before considered new walk   
    #adj accelth to 1.4 from 1.5 after first data drop from radbout 
    accelTh = 1.4 # threshold of aceleration to determine whether walking related 
    accelThhi = 7 # upper threshold of aceleration to determine whether walking related 
    walkGapMaxs = 2 # max  gap in seconds  between steps before considered new walk   

    maxBreachLen = 10 # can't be accelerating vertically for very long ! 
    maxBreachLens = 0.25 # can't be accelerating vertically for very long ! 
    minBreachLen = 3 # potential breach becomes a breach here 

    minBreachLens = 0.05 # potential breach becomes a breach here 
    step = {} # dict to record time of centre of breach for each step in a walk 


class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    

def isendofBreach(i,y):
    ''' are we at the end of the breach ? need > blipsize readings < threshold
    if signal at i-blipsize:i < thresh endof = True 
    also returns how long ago breach ended
    
    '''
    global accelTh,blipSize
    if  all([x < accelTh for x in y[i-blipSize:i+1]]):
       
        return True, i-blipSize
    else: 
         print(i,' blip detected')
         return False, 0
     
def getWalkAcc(walklog,data,walkno,samplerate):
    '''for a given walk no in a wlk log , get accelerometer data '''
    twstart = walklog[walkno][0]
    twend = walklog[walkno][1]
    twstart *= samplerate
    twend *= samplerate
    twstart = int(twstart)
    twend = int(twend)
    axyz = pd.DataFrame(data[twstart:twend])
    
    return axyz

def recordEndWalk(walkCount,stepFirst,stepLast,stepCount,avact,epoch,stepstats):
    # startwalk in epoch format 
    global walklog
    wst = timedelta(seconds=stepFirst)
    wst = pd.to_datetime(epoch[0]) + wst
    walklog[walkCount] = [stepFirst,stepLast,stepCount,avact,wst,stepstats]
    return
    
def getWalks(y,actclasses,smawindow,epoch,samplerate=50,accelthov=1.4,statslen=25):
    ''' get walks from input vertical acceleration data 
    
    activity type - the activity type can be 0, no activity (sit),
        1 medium(walk), 2 high(run) , 3 very high (potential fall). 
        They are measured over a window, usually of one second.
        rather than exclude walks with high intensity , will mark as
        such and analyse later 
        input: 
            statslen - we will generate stats every statslen steps for use by 
            gait routine 
        
        output: walklog - list of walks, strt (secnds offset ), end (seconds offset), stes, activity levels 
                timestamps (tod) . 
    
    '''
    # x counts seconds into file, i counts positions 
    #y = vacc
    x = np.arange(0,len(y))/samplerate # smm is 50-51 Hz 
    x = x-x.min()
   
    initValues() # warning - lots of globals here walk detection 
    global potentialBreach,breachState,breachStart,walkState,walkCount, \
     walklog,       \
     stepCount,     \
     walkCount,     \
     blipSize,      \
     stepGapMin,    \
     walkGapMin,    \
     accelTh,       \
     accelThhi,       \
     maxBreachLen,  \
     minBreachLen,  \
     step
    accelTh = accelthov
    # all times in seconds - converted by sectopoints
    global walkGapMax
    global stepGapMins, walkGapMaxs, maxBreachLens, minBreachLens
    stepstats=[] # record stats every statlen steps in here 
    stepstatcounter = 0 # count the number of 25 step segments 
    stepdur = {}
    stepstatidx = {} # indices for 25 step stats segments
    # y = simplea.copy() # sample data 
    def getAvAct(samplestart,sampleend,actclasses,smawindow):
        ''' calculate average activity level period '''  
        actstart = int(round(samplestart/smawindow,None))
        actend = int(round(sampleend/smawindow,None))
        print(actstart,actend+1,type(actstart),type(actend))
        return np.mean(actclasses[actstart:actend+1])
 
    #TODO incorporate upper breach limit here .... if exceeded, counts as end of breach 
    # as we've moved into someother activity 
    for i in np.arange(len(y)):
        if y[i] > accelTh:
            # if we're already in breach: fail if > x timesteps have passed (FAIL)
            if breachState:
                if x[i]-breachStart > maxBreachLen:
                    # note in my walking a breach lasted about 1/6 sec 
                    raise MyError(i,'gone mad - more than maxbreachlength readings gt threshold')
                
            else:
               # if y[i] > accelThhi: continue # exceeded uper threshold 
                print(i,'breachstart, time ',x[i])
                breachState = 1
                breachStart = i # may need to account for blips here t
                breachStarts = x[i] # may need to account for blips here t
                # but should be taken care of by minbreach length 
                
        # !! Fallen below threshold !!                                  
        elif breachState: # y[i] < thresh - possible end of breachstate
            endofBreach, breachEnd = isendofBreach(i,y)
            breachEnds = x[breachEnd]
            if endofBreach:
                print(i,'end of breach detected')
         
                # if it's not long enough, chuck it away 
                if breachEnds - breachStarts < minBreachLens:
                    # ok we have a min breach - we will ignore 
                    # ... for now 
                    print(i,'mini-breach - chuck it. time ',x[i])
                    breachState = 0
                    continue # get next point 
                # record midpont of step        # calculate breachCentre
                step[stepCount] = breachStarts + (breachEnds-breachStarts)/2 
                #record duration of step 
                
                breachState = 0 
                # calc previous step length (if there was one)
                # we only know walk has begun at second step ... or end of first step 
                # 1st step is complete at end of minstepgap  time after end of breach 
                if walkState: # in a walk already - add new step 
                
                    stepLen = step[stepCount]-step[stepCount-1]
                    stepdur[stepCount] = stepLen
                    print(i,' steplen: ',stepLen)
                    # if gap is less than minstep gap, something is wrong 
                    if stepLen < stepGapMins:
                        print(i,'stepLen ', stepLen,' < stepGapMin',stepGapMins,' - not expected')
                    elif stepLen < walkGapMaxs: # gap between steps OK  - part of current walk
                        stepCount += 1
                        if stepCount%statslen == 0: # if multiple of statslen, generate stats
                        # start of segment, end, sd of durations
                            if stepstatcounter == 0:
                                stepstatidx[0] = breachStart
                            stepstatcounter += 1 
                            stepstatidx[stepstatcounter] = i # sample at which this segment ends 
                            stepdurs  = [v for k,v in stepdur.items()   \
                            if ((k>=stepCount-statslen) & (k<stepCount))]
                            stepstats.append([stepstatidx[stepstatcounter-1],i,
                                             np.asarray(stepdurs).std()])
                             
                           
             
                else: # start of a new walk 
                    print(i,'start of walk')
                    walkState = 1
                    walkCount += 1
                    stepCount += 1
        elif walkState: # not breach state - continuation of below threshold 
            
            # if time since last breach greater than maxwalk threshold , delcare end of walk 
            if x[i]-step[stepCount-1] > walkGapMaxs: # max gap allowed in a walk 
                        print(i,'end of walk')
                        # get mean activity level for walk 
                        avact = getAvAct(step[0]*samplerate,step[stepCount-1]*samplerate,actclasses,smawindow)
                        recordEndWalk(walkCount,step[0],step[stepCount-1],stepCount-1,avact,epoch,stepstats)
                        walkState = 0 
                        stepCount = 0 
            else:
                pass
               # print(i, 'not at walkGapMaxs ',x[i]-step[stepCount-1])
    #if we have reached the end and we'rein the middle of a walk, better record it
    if walkState:
        avact = getAvAct(step[0]*samplerate,step[stepCount-1]*samplerate,actclasses,smawindow)
        recordEndWalk(walkCount,step[0],step[stepCount-1],stepCount-1,avact,epoch,stepstats)
        walkState = 0 
        stepCount = 0 
        
    return walklog
#%% postural transition functions 
def getVertVelocity(vacc,plot=0):
    ''' requires vertical body acceleration, output from getVectors
    '''
    # filter vertical acceleration to remove drift. this does not add a lag 
    Vaccf =  butter_bandpass_filter(vacc,0.15,15,50,2)
    #Vaccf =  butter_bandpass_filter(arb,0.15,15,50,2)
    x=np.arange(0,len(Vaccf),1)/50 # using dx instead 
    Vv = cumtrapz(Vaccf,dx=1/50,axis=0,initial=0)
    
   # Vvr = cumtrapz(acc,x=np.arange(0,len(acc),1)/50)
    if plot:
        plt.plot(Vv,label='vert velocity')
        #plt.plot(Vaccf[600:]/9.8,label='vert accel') # scales nicely measured in g 
        plt.legend()
        plt.plot()
    return Vv

def getSintheta(acc,sma):
    ''' 
    input : acceleromemter inc. gravity  data N x 3 
    loop through sma (1 second windows) < 2m/s (curone 2010)
    note filtering parameters set here 
     '''
    #for now, use first 600 to gather static a 
    # will use SMA to do this automatically
    # note: we are using the smoothed accelerometer data, no filters 
    # applied 
    # automating quasi static 
    # choosing min length 5 second of < 1.5m/s for definition of quasi
    # static event 
    # set theta 0 , keep at zero until next 5 sec inactive segment 
    iasegs = getInactiveSegments(sma,5,1.5)
    #for each inactive segment, set theta = zero 
    # 2. calculate a
    # for subsequeny non active segment calc theta using vector algebra 
    # contine until end 
    sintheta = []
   
     # need to add in the value for when the sequence starts off active
    # segment counts start at 1 
    if iasegs[0][0] > 1:
        #fill with zeros for first x segs 
        sintheta = np.zeros((iasegs[0][0]-1)*50)
    for i,iaseg in enumerate(iasegs):
        startseg = (iaseg[0]-1)*50 # iaseg numbered from 1 
        endseg = (iaseg[1]-1)*50 # deliberately leaving out last 1s segment 
        samples = (endseg - startseg) # no. samples inseg given segs are 1 second each 
        sintheta = np.concatenate([sintheta,np.zeros(samples)])
        Va =  np.mean(acc[startseg:endseg],axis=0)
        # get start of next seg or end of sma 
        if i+1 == len(iasegs):
            print('reached end of iasegs - this s expected')
            endactive = len(acc)
        else: #ends at the beginning of next inactive segment 
            endactive = (iasegs[i+1][0]-1)*50
        # calculate theta from endseg to endactive 
        Vb = acc[endseg:endactive]  # no need for this ? 
        Vb = Vb.transpose()
        VadotVb = np.dot(Va,Vb)
        magb = np.square(Vb)
        magb = magb.sum(axis=0)
        magb = np.sqrt(magb)
        maga = np.sqrt(np.sum(np.square(Va)))
        magba = magb*maga
        costhetaseg = VadotVb/magba
        sinthetaseg = np.sqrt(1-np.square(costhetaseg))
        sintheta = np.concatenate([sintheta,sinthetaseg])
        
    return sintheta

def synchWalks(walklog,epoch):
    ''' add synched time to walk logs
    input : wlaklog and epoch indicating starttime ''' 
    
    return swalklog

def movingAv(vec,window_width):
    ''' calculate moving average
        # thanks to https://stackoverflow.com/users/2370124/roman-kh for this super-simple solution 
    '''
    cumsum_vec = np.cumsum(vec) 
    #window_width=12
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    pad1 = window_width
    pad2 = 0
   #pad1 = window_width//2
   # pad2 = window_width-pad1
    return np.pad(ma_vec,(pad1,pad2),'edge')

def getPT(walklog,sintheta,rate,angleth=0.35,excludewalk=0):
    # count potential PTs as start sintheta > angleth, end < angleth
    # perioods identified as walking can be excluded from the data or not 
    PTstate = 0 # whether we have detected a PT state 
    #create boolean from walk data - 1 for inwalk, 0 for not in walk 
   
    walkbool = np.zeros(len(sintheta),dtype=int)

    walks = [v  for k,v in walklog.items() if v[2] > 2]
   
    for walk in walks:
        wkstart = int(round(walk[0]*rate))
        wkend = int(round(walk[1]*rate))
        wklen = wkend-wkstart
        walkbool[wkstart:wkend] = np.ones(wklen)
       
    # remove offset 
    #walkbool = walkbool[offset:]
     
    PTlog={}  
    PTcount=0
    PTstate=0
    PTmax = 2*rate # no transition > 2 sec 
   # sinfilt = medfilt(sintheta,7)

    sinfilt = movingAv(sintheta,12)
    # TODO : need to allow for changes in sintheta, not absolute threholds. 
    # so if sin theta falls to 50% of the filtered peak, and then goes back up , 
    # even if the valley is < threshold, counts as a new peak ! 
    # may be better to look at itegral of sintheta to allow for drops 
    # below the peak 
    i=-1
    for theta in sinfilt:
        i+=1
        if theta > angleth:
            if PTstate:
                # TODO if i-PTstart > PTmax:
                    #end walk
                    #increment i until theta < angleth
                continue
            else:
                if (excludewalk) & (walkbool[i]): # ignore if part of walk 
                    print(i,' part of walk')
                else:
                    PTstate=1
                    PTstart = i
        else:
            if PTstate:
                PTstate = 0 
                PTlog[PTcount] = [PTstart,i,i-PTstart]
                PTcount += 1
    #TODO amalgamate peaks where gap < x milliseconds and average from start of one to end of next is > thresh 
                
    return PTlog

def getPTdetail(vacc,sintheta,PTlog,smph=0.35,smpd=40):
    ''' get the type of postural transition : 
        0 - none: not a PT - don't know what it is ! 
        1- sist
        2 - stsi
        3 - bumpy sist (double and or twist)
        4 - bumpy stsi
        and add duration  
        input: VERTICAL  accelerometer data as array, log of start end  of PTS
        (PTlog) in sample times,
        smph: minimum peak height for sintheta
        smpd: minimum peak distance for sintheta - in samples '''
    minlength = 20 # samples = /50 seconds * 50Hz - minimum duration of PT      
    from BMC.functions import detect_peaks
    Vv = getVertVelocity(vacc)
    sinfilt = medfilt(sintheta,25) # this may be too high 
    sinfilt = movingAv(sintheta,12)
   #plt.plot(Vv,label='Vv')
    PTdetail = {}
   # len(sintheta)
    for i,pt in PTlog.items():
        start,end,dur = pt 
        #start,end,dur = PTlog[33] # testing 
        if dur < minlength:
            continue
        showgr = False
        #plt.plot(sinfilt,label='sinfilt')
       # addGrid()
        ptsintheta = sinfilt[start:end]
        plt.plot(ptsintheta)
        # now we have one PT. get peaks for sintheta, then peaks for arb 
        sinthpeaks = detect_peaks.detect_peaks(ptsintheta,mph=smph,
                                               mpd =smpd,
                                               edge='both',
                                               kpsh=True,
                                               show = showgr)
       
        # if 2 take the mean, else the median ?
        # TODO - if there are multiple peaks x distance apart, raise hte threshold to 
        # separate them 
        if len(sinthpeaks) == 0:
            print(i,'no sin theta peaks found, rejecting PT. mph,mpd =',smph,smpd)
            continue
        try:
            sinthpeak = np.mean(sinthpeaks)
            sinthy = ptsintheta[int(sinthpeak)]
        except ValueError:
            print('sinthpeaks: ',sinthpeak, 'PTlog: ',i)
            continue
        #plt.plot(ptvv)
        #addGrid()
         #get 1 secs before and after peak 
        sinpeak = int(start+sinthpeak)
        offset = sinpeak-50 # sample no. of start of data 

        ptvv = Vv[offset:sinpeak+50]
        Vvpeaks = detect_peaks.detect_peaks(ptvv,mph=0.25,
                                            mpd = 40,
                                            show = showgr)
        if len(Vvpeaks) == 0:
            print(i,' no Vertical velocity peaks found, rejecting PT. mph=',)
            continue
        
        Vvvalleys = detect_peaks.detect_peaks(ptvv,mph=None,mpd=20,
                                              valley=True,
                                              show = showgr)
        if len(Vvvalleys)==0:
            print(i,' no valleys found, rejecting PT')
            continue
        # we want the valley before the peak 
        Vpeak = np.mean(Vvpeaks) # assuming max of 2? may need the one closest to sintheta peak 
        # find smallest negative value to get closest valley before peak 
        Vvv = Vvvalleys - Vpeak # this is a float , so we can:
        idx = (Vvv**-1).argmin()
        Vvalley = Vvvalleys[idx]
        # how far apart are peak and valley ... min 80, max 220 (by eye)
        #TODO reject if samples too far apart 
        #note that the velocity zero can drift around , so if the valley is above zero, 
        # adjust to -0.1 and make same adjust to vpeak 
        if ptvv[Vvalley] >= 0:
            valley = -0.1
            diff = ptvv[Vvalley] - valley
            peak = ptvv[int(Vpeak)] - diff 
        else:
            valley = ptvv[Vvalley]
            peak = ptvv[int(Vpeak)]
            
            
        ratio = abs(peak/valley)
        if ratio > 1.9: result = 1 #sist ratio arrived at empirically 
        #TODO derive the ratio from the data 
        else: result = 2 # stsi
        possintheta = [sinpeak,sinthy]
        pospeak = [Vpeak+offset,peak]
        posvalley = [Vvalley+offset,valley]

        #returnval  = [possinthteta,pospeak,posvalley]
        PTdetail[i] = [result,possintheta,pospeak,posvalley,ratio] # will add the duration later 
        
             
    return PTdetail



#%% Gait analysis

def getGaitfeature(axyz,gravity,samplerate,row,plot_test=False):
    ''' get gait feature for one walk
    output: rc - 1 is valid, anything else not 
            features df from extractGait        
    '''
    axyz /= gravity # mhealthx expects values in g not ms-2
           # plt.plot(vacc18,label='vert acc')
           #addGrid() 
           # row to prefix to features generated by extract gait (it's the way the mhealth code works) 
    rc, features =  extractGait(axyz,samplerate,row,plot_test)
    return rc, features 

def getGaitfeatures(walklog,arb,walkno='all',samplerate=50,gravity=9.81,plot_test=False): 
    ''' input : 
        log of start and end of walks 
        accelerometer array, body vectors only 
        
        '''
        #walkno=47
        #walklog = walks.copy()
   # del features
    import pandas as pd
 
    from mhealthx.extract import run_pyGait
    import mhealthx.extractors.pyGait
    import mhealthx.signals
    from mhealthx.extractors.pyGait import project_walk_direction_preheel, walk_direction_preheel
    # imp.reload(mhealthx.extractors.pyGait)
    
    featuresdf = pd.DataFrame()
    walkerrs = []
    
    for walk,detail in walklog.items():
        if (walkno != 'all') & (walkno != walk): continue
       
        
        twstart, twend, steps, act,dtime, stepstats  = detail
  #
  
        
        #  walk=49; twstart, twend, steps, act, dtime = walks[walk]
        # note, for ADL we should perhaps be looking at 20 steps minimum  - see
        # Published online 2012 Feb 8. doi:  10.1186/1743-0003-9-11
        
        if (steps > 7) & ((act > 0.5) & (act < 1.8)):  # walklength(no. steps)  > 2 
            print('start walk ',walk,steps,act)
            # remove first and last steps 
            # step duration = walk duration/no. steps 
            stepduration = (twend-twstart)/steps
            twstart = twstart + stepduration
            twend = twend - stepduration
            twstart *= samplerate
            twend *= samplerate
           
            # the mhealthx gait extraction routines cannot handle any noise, 
            # so I have decided not to attempt to give it the first and last steps - 
            # in terms of gait, it is good to ignore the first steps anyway
            # and last steps anyway 
            # as per the 6 or 10 mw tests 
           # twstart = twstart-(int(round(samplerate/4,0)))  # walk identification starts from  centre of step
          #  twend = twend + (int(round(samplerate/4,0)))  # gait function needs start of step + lead in/out
            twstart = int(twstart)
            twend = int(twend)
            print('twstart,twend',twstart,twend)

            axyz = arb[twstart:twend,:].copy()
           # pv = vacc[twstart:twend]
       #     pv= pv/gravity
         #   plt.plot(axyz)
        #   vacc18 = vacc[twstart:twend]
            # axyz=vex.copy()
            row = pd.Series({'walkid':walk, 'start':twstart, 'Method1_stepcount':steps})           
            rc, features =  getGaitfeature(axyz,gravity,samplerate,row,plot_test=plot_test)
            if rc != 1:
                print('bad rc from walk:',walk)
                walkerrs.append(walk)
            else:
                print(features)
                featuresdf = pd.concat([featuresdf,features])
                
        else: print('walk ',walk,' rejected',steps,act)
        
    return featuresdf, walkerrs

   



def extractGait(axyz,samplerate,row,plot_test):
    # input_file = cfg.WRITE_PATH + '/testfiles/' + 'accel_walking_return-test-1.csv'
   # basePath = "C:\\Users\\dad\\.synapseCache\\428\\5389428\\428\\5389428\\"

    #tf2 = "accel_walking_outbound.json.items-502d0501-43ee-4462-a7e4-012e216576922296422479994601681.tmp"
   # input_file = basePath + tf2
  
  
#    device_motion = False
#    t, axyz, gxyz, uxyz, rxyz, samplerate, duration = read_accel_json(input_file, start, device_motion)
#    if samplerate == 0: return -1,'File of too short duration'
    # axyz is body motion only see  Mpower docs .... or is it -need to get the gravity component, so doubt it
    ax, ay, az = axyz.transpose()
    #ax, ay, az = vex[5:].transpose()
    
    duration = len(ax)/samplerate
    t=list(np.arange(0,len(ax)/samplerate,1/samplerate))
    # at times the clock seems to vary, so sample rate may not apply to this segment
    # make it same length as segment ax 
    t = t[:len(ax)]
    
    #type(axyz)
    stride_fraction = 1.5/8.0 # MT changed from 1/8]
    threshold_heel_strike = 0.25
    gaitthreshold = 0.25
    order = 4
    cutoff = max([1, int(samplerate/10)])
    distance = None
    #row = pd.Series({'a':[1], 'b':[2], 'c':[3]})


    save_rows = True # writes to csv rater than synapse table 
    #walk_direction_preheel
    #def walk_direction_preheel(ax, ay, az, t, samplerate, 
    #                          stride_fraction=1.0/8.0, threshold=0.5,
    #                         order=4, cutoff=5, plot_test=False)
   
    

   # directions = walk_direction_preheel(ax, ay, az, t, samplerate,\
   #                                             stride_fraction, threshold0, order, cutoff)
 
    # imp.reload(mhealthx.extract)
    # imp.reload(mhealthx.signals)

    # imp.reload(mhealthx.extractors.pyGait)
    from mhealthx.extractors import pyGait
    from mhealthx.extract import run_pyGait
    from  mhealthx.signals import autocorrelate
    file_path='.' # not using 
    table_stem=''
    print('samplerate:',samplerate)
    try:
        px, py, pz = project_walk_direction_preheel(ax, ay, az, t, samplerate,
                                                stride_fraction, threshold_heel_strike, 
                                                order, cutoff,plot_test=plot_test)
        
        feature_row = run_pyGait(px, t, samplerate,
                                            duration, gaitthreshold, order,
                                            cutoff, distance, row, file_path,
                                            table_stem, save_rows,plot_test=plot_test)
      #  plt.plot(np.array([px,py,pz]).T,label=['x','y','z'])
       # plt.legend()
    except (ValueError, IndexError) as e:
        return -1,e
    else:
        return 1,feature_row
    

                




        