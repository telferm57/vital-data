# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 18:05:31 2017

@author: telferm
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:12:47 2017

@author: telferm
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.signal import medfilt
from scipy.integrate import cumtrapz 
import imp 
from math import floor
import pandas as pd
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
    from vitalsensor import butter_lowpass
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=5):
    from vitalsensor import butter

    nyq = 0.5 * fs
    cutf = cutoff / nyq
    b, a = butter(order, cutf, btype='low',analog=False)
    return b, a

def lpFilter(x,y,z,cutoff,samplerate,order=4):
    from vitalsensor import butter_lowpass_filter
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
    ar2 = np.square(art)
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
    return ar.max()-ar.min()

def testOffBody(data):
   
    ''' given input accelerometer data, return whether device was off body or not
     #      https://doi.org/10.1371/journal.pone.0022922 Van Hees 2011
     
     2 axes have a std deviation < 0.1; OR 2 axes have a range < 0.8 
     
    " A block was classified as non-wear time if the standard deviation was
     less than 3.0 mg (1 mg = 0.00981 m·s−2) for at least two out of the three
     axes or if the value range, for at least two out of three axes, was less 
     than 50 mg"
     
     3mg = 0.03 m/s 50 mg = 0.1 m/s 
    '''
    def getRange(ar):
        return ar.max()-ar.min()
    
    rangethr = 0.8; sdthr=0.1
    res = 0 # it IS on body unless ...
    axsdtot = 0; axrangetot = 0 
    for col in data.columns:
       
        sdaxis = np.std(data[col].values)
        if sdaxis < sdthr: 
            #print(col,' sd gt 0.03 :',sdaxis)    
            axsdtot += 1
        print(col, sdaxis, getRange(data[col].values))           
        if getRange(data[col].values) < rangethr: 
            #print(col,' sd gt 0.03 :',sdaxis)    
            axrangetot += 1
    if (axsdtot > 1) |  (axrangetot > 1) : res = 1         
    print(col,' no.  axes exceeded 1 :',axsdtot,axrangetot)    
    return res

def getOffBody(ws,data,rate):
    ''' for a given window size, return a dict of windows indicating whether device off or 
    on body '''
#    data=tt.copy()
#    ws = 2 # minutes
#    rate=50
    offbody = {}
    samplesize = ws*rate*60 # 6000 per 2  minutes 
    for sample in np.arange(floor(len(data)/samplesize)):
        sstart = sample*samplesize
        send = sstart+samplesize
        offbody[sample] = testOffBody(data[sstart:send])
        
    return offbody  

def dissectData(data,onbodywindow=2,samplerate=50):
    ''' for a given datastream, identify on-body segments and pickle them 
    input: pandas dataframe of data with accelerometer in columns named x,y,z
            window size in minutes 
    
    returns number of periods dict '''
    
    tt= data.copy()
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
    np.std(tt.z.values)
    tt.describe()
    tt.drop(['x','y','z'],axis=1,inplace=True)
    tt.plot()
    offBodySegments = getOffBody(onbodywindow,tt,samplerate) # x minute mindow, 50 samples/second 
   #plt.plot(np.asarray(list(offBodySegments.values())))
    #plt.show()
    #only interested in > 10 minutes of wear time, so 5 consecutive zeros  : 
    ar = np.asarray(list(offBodySegments.values()))
    # get periods of segments from ' off body' dict 
    started=0
    runsize = 6
    periodcounter = 0 
    periods = {}
    for i in np.arange(len(ar)):
        
        if ar[i:i+runsize].sum() == 0:
           
            if started:
                pass
            else:
                startrun = i
                print('start at ',i)
                started = 1
                periodcounter += 1
        else:
            if started: #stop
                started = 0 
                endrun = i+runsize-2
                print('end: ',endrun)
                periods[periodcounter] = [startrun,endrun]
    
    #create seperate files from periods and pickle them . periods are stored as 
    # sample numbers 1-n , so is dependent on sample size set for the run 
    ws=2 # minutes
    samplesize = ws*samplerate*60
  
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
                
    fn='userid-date-device'            
    picklePeriods(data, periods, samplesize,fn)
   # import pdb
    #pdb.set_trace()
    
    return periods, offBodySegments


def getGravity(data,offBodyPeriods,onbodywindow):
    ''' gravity reading will vary for each device.. this routine will calculate
    what the particular device is measuing as 1g. that value will be used 
    throughout to calculate thresholds. off body periods can be used for tis 
    
    '''
    gravlist = []
    for key, val in offBodyPeriods.items():
        if val == 1:
            start = key*samplerate*onbodywindow*60
            end = start + samplerate*onbodywindow*60 -1
            ar = smm[start:end].values
          
           
            gmag = magnitude(ar.transpose())
            gest = np.mean(gmag)
            gravlist.append(gest)
            #print('period ', key, gest)
            
    
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
    from vitalsensor import lpFilter
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
    newdf = data.iloc[7:,:]
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
     
     output : list of inactive segments greater than minlength log.
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
     walklog,       \
     stepCount,     \
     walkCount,     \
     blipSize,      \
     stepGapMin,    \
     walkGapMin,    \
     accelTh,       \
     maxBreachLen,  \
     minBreachLen,  \
     step
    # all times in seconds - converted by sectopoints
    global walkGapMax
    global stepGapMins, walkGapMaxs, maxBreachLens, minBreachLens
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
    
def recordEndWalk(walkCount,stepFirst,stepLast,stepCount,avact):
    walklog[walkCount] = [stepFirst,stepLast,stepCount,avact]
    return

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


    
def getWalks(y,actclasses,smawindow,samplerate=50):
    ''' get walks from input vertical acceleration data 
    
    activity type - the activity type can be 0, no activity (sit),
        1 medium(walk), 2 high(run) , 3 very high (potential fall). They are measured over a window, 
        usually of one second . rather than exclude walks with high intensity , will mark as such 
        and analyse later 
    
    '''
    # x counts seconds into file, i counts positions 
    x = np.arange(0,len(y))/samplerate # smm is 50 Hz 
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
     maxBreachLen,  \
     minBreachLen,  \
     step
    # all times in seconds - converted by sectopoints
    global walkGapMax
    global stepGapMins, walkGapMaxs, maxBreachLens, minBreachLens
    
    # y = simplea.copy() # sample data 
    def getAvAct(samplestart,sampleend,actclasses,smawindow):
        ''' calculate average activity level period '''  
        actstart = int(round(samplestart/smawindow,None))
        actend = int(round(sampleend/smawindow,None))
        print(actstart,actend+1,type(actstart),type(actend))
        return np.mean(actclasses[actstart:actend+1])
 
    
    for i in np.arange(len(y)):
        if y[i] > accelTh:
            # if we're already in breach: fail if > x timesteps have passed (FAIL)
            if breachState:
                if x[i]-breachStart > maxBreachLen:
                    # note in my walking a breach lasted about 1/6 sec 
                    raise MyError(i,'gone mad - more than maxbreachlength readings gt threshold')
                
            else:
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
                # calculate breachCentre
                # if it's not long enough, chuck it away 
                if breachEnds - breachStarts < minBreachLens:
                    # ok we have a min breach - we will ignore 
                    # ... for now 
                    print(i,'mini-breach - chuck it. time ',x[i])
                    breachState = 0
                    continue # get next point 
                # record midpont of step 
                step[stepCount] = breachStarts + (breachEnds-breachStarts)/2 #record length of step (s)
                breachState = 0 
                # calc previous step length (if there was one)
                # we only know walk has begun at second step ... or end of first step 
                # 1st step is complete at end of minstepgap  time after end of breach 
                if walkState: # in a walk 
                    stepLen = step[stepCount]-step[stepCount-1]
                    print(i,' steplen: ',stepLen)
                    # if gap is less than minstep gap, something is wrong 
                    if stepLen < stepGapMins:
                        print(i,'stepLen ', stepLen,' < stepGapMin',stepGapMins,' - not expected')
                    elif stepLen < walkGapMaxs: # gap between steps OK  - part of current walk
                        stepCount += 1
                        
             
                else: # start of a new walk 
                    print(i,'start of walk')
                    walkState = 1
                    walkCount += 1
                    stepCount += 1
        elif walkState: # not breach state - continuation of below threshold 
            
            # if time since last breach greater than maxwalk threhold , delcare end of walk 
            if x[i]-step[stepCount-1] > walkGapMaxs: # max gap allowed in a walk 
                        print(i,'end of walk')
                        # get mean activity level for walk 
                        avact = getAvAct(step[0]*samplerate,step[stepCount-1]*samplerate,actclasses,smawindow)
                        recordEndWalk(walkCount,step[0],step[stepCount-1],stepCount-1,avact)
                        walkState = 0 
                        stepCount = 0 
            else:
                print(i, 'not at walkGapMaxs ',x[i]-step[stepCount-1])
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
    acc=ars.copy() # temporary for testing
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
            print('end of iasegs detected')
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

def getPT(walklog,sintheta,rate,angleth = 0.35):
    # count potential PTs as start sintheta > angleth, end < angleth
    # can as a boolean vector with walk as the same - later maybe 
    PTstate = 0 # whether we have detected a PT state 
    #create boolean from walk data - 1 for inwalk, 0 for not in walk 
    offset = 0 
    walkbool = np.zeros((endperiod-startperiod)+offset,dtype=int)
     # get min max for walk speeds to get rgb limits 
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
    sinfilt = medfilt(sintheta,7)
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
                if walkbool[i]: # ignore if part of walk 
                    print(i,' part of walk')
                else:
                    PTstate=1
                    PTstart = i
        else:
            if PTstate:
                PTstate = 0 
                PTlog[PTcount] = [PTstart,i,i-PTstart]
                PTcount += 1
                
    return PTlog

def getPTdetail(vacc,sintheta,PTlog):
    ''' get the type of postural transition : 
        0 - none: not a PT - don't know what it is ! 
        1- sist
        2 - stsi
        3 - bumpy sist (double and or twist)
        4 - bumpy stsi
        and add duration  
        input: accelerometer data as array, log of start end in sample times '''
    minlength = 20 # samples = /50 seconds * 50Hz - minimum duration of PT      
    from BMC.functions import detect_peaks
    Vv = getVertVelocity(vacc)
    sinfilt = medfilt(sintheta,25) # this may be too high 
    plt.plot(Vv,label='Vv')
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
        # now we have one PT. get peaks for sintheta, then peaks for arb 
        sinthpeaks = detect_peaks.detect_peaks(ptsintheta,mph=0.4,
                                               mpd = 40,
                                               edge='both',
                                               kpsh=True,
                                               show = showgr)
       
        # if 2 take the mean, else the median ? 
        sinthpeak = np.mean(sinthpeaks)
        sinthy = ptsintheta[int(sinthpeak)]
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
            continue
        
        Vvvalleys = detect_peaks.detect_peaks(ptvv,mph=None,mpd=20,
                                              valley=True,
                                              show = showgr)
        if len(Vvvalleys)==0:
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
#%% plotting helpers
def addGrid():
    ax = plt.gca()
    ax.minorticks_on()
    ax.grid(which='both')
    ax.legend()

def plotWalkTimes(ax, walklog,rate=50, param_dict={}):
    ''' given a walklog object, that records all walks over a predfined sensor
    epoch, and a set of axes, draw the blocks of time on the axes 
    corresponding to walks > 2 steps. walks < 2 steps may be something else 
    
    input: walklog - dict of walks that have start& end  time (seconds) and
    number of steps per walk 
    '''
  
    # get min max for walk speeds to get rgb limits 
    ll = [(v[1]-v[0])/v[2]  for k,v in walklog.items() if v[2] > 2]
    # scale min-max for g=rgb 0.25 - 0.75
    vrange = max(ll) - min(ll)
    vmin = min(ll)
    ll2 = ((ll-vmin)/vrange)*0.5+0.25
    
    
    for walk,detail in walklog.items():
        if detail[2]> 2:  # walklength(no. steps)  > 2 
            
            x = detail[0]*rate
            y= -10
            w= (detail[1]-detail[0])*rate
            h=20
            #x = 20000
            y = -10
           # w=20000
            steptime= w/(rate*detail[2])
            rgbrval = ((steptime-vmin)/vrange)*0.8+0.2
            ax.add_patch(patches.Rectangle((x,y),w,h,
                                           fc=(rgbrval,0,0)))
            
#            ax.text(x, y, 'center top',
#                    horizontalalignment='center',
#                    verticalalignment='top',
#                    transform=ax.transAxes)
            print(x,w)
            ax.text(x,10,detail[2])
            ax.text(x,9.5,'{0:1.2f}'.format(steptime),color='lightblue',
                    rotation=90)        
    out = ax.plot()
    return out

def plotPTTimes(ax, PTlog,rate=50,offset=0,param_dict={}):
    '''  plot pt times on axes  '''
    h=20
    # PTs last between 0.5 and 4 seconds .. so throw away anything 
    # 1/2 the rate and 4 * rate- although are PTs contigous ? they are in FTSTS 
    for PT in (y for y in PTlog.values() if (y[2] >= 20)):
        x,x2 = PT[0],PT[1]       
        w= x2-x +1
     
        xoff = x  + offset # 600 used to calc quasi-static
        print(x,w,xoff)
        y=-10
        ax.add_patch(patches.Rectangle((xoff,y),w,h,
                                           fc='b',alpha=0.2))
        
    out = ax.plot
    return out

def plotPTdetail(ax,PTdetail,rate=50,param_dict={}):
    ''' plot PT details on axes - intended to overlay the Vacc and sintheta 
    plots 
    
    PT detail consists of the classification of the PT, 3 x,y 
    co-ordinates of the sintheta peak and the vert velocity trough and peak 
    and the ratio of peak to trough. The ratio is used to classify the 
    transition '''
    # get first set of points
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
                

    
    out = ax.plot

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


#%% Pressure Data
def getPressureData(fn,nrows=1e4,skiprows=0):
    import pandas as pd
    path='C:/Users/telferm/python/projects/walkingsignals/walkingsignals/data/smm/pilot/' 
    smm = pd.read_csv(path + fn,header=None,skiprows=skiprows,nrows=nrows)
    return smm 

#%% Gait analysis

def getGaitfeatures(walklog,ar): 
    from mhealthx.extract import run_pyGait
    import mhealthx.extractors.pyGait
    import mhealthx.signals
    from mhealthx.extractors.pyGait import project_walk_direction_preheel, walk_direction_preheel
    # imp.reload(mhealthx.extract)
    
    samplerate = 50
    
    featuresdf = pd.DataFrame()
    
    for walk,detail in walklog.items():
        
        twstart, twend, steps, act  = detail
        #twstart, twend, steps, act = walklog[15]
        if (steps > 10) & ((act > 0.5) & (act < 1.8)):  # walklength(no. steps)  > 2 
            print('start walk ',walk,steps,act)
            twstart *= samplerate
            twend *= samplerate
            twstart = int(twstart)
            twend = int(twend)
            axyz = ar[twstart:twend]
        #   vacc18 = vacc[twstart:twend]
            axyz /= 9.81
           # plt.plot(vacc18,label='vert acc')
           #addGrid() 
           # row to prefix to features generated by extract gait (it's the way the mhealth code works) 
            row = pd.Series({'walkid':walk, 'start':twstart, 'Method1_stepcount':steps})
            features =  extractGait(axyz,samplerate,row)
            featuresdf = pd.concat([featuresdf,features[1]])
        else: print('walk ',walk,' rejected')
        
    return featuresdf

   



def extractGait(axyz,samplerate,row):
    # input_file = cfg.WRITE_PATH + '/testfiles/' + 'accel_walking_return-test-1.csv'
   # basePath = "C:\\Users\\dad\\.synapseCache\\428\\5389428\\428\\5389428\\"

    #tf2 = "accel_walking_outbound.json.items-502d0501-43ee-4462-a7e4-012e216576922296422479994601681.tmp"
   # input_file = basePath + tf2
    start = 0
  
#    device_motion = False
#    t, axyz, gxyz, uxyz, rxyz, samplerate, duration = read_accel_json(input_file, start, device_motion)
#    if samplerate == 0: return -1,'File of too short duration'
    # axyz is body motion only see  Mpower docs .... or is it -need to get the gravity component, so doubt it
    ax, ay, az = axyz.transpose()
    
    duration = len(ax)/samplerate
    t=list(np.arange(0,len(ax)*0.02,0.02)[:-1])
    
    type(axyz)
    stride_fraction = 1.0/8.0
    threshold0 = 0.5
    gaitthreshold = 0.4
    order = 4
    cutoff = max([1, samplerate/10])
    distance = None
    #row = pd.Series({'a':[1], 'b':[2], 'c':[3]})

    file_path = 'test1.csv'
    table_stem = './data'
    import os
    os.path.isdir(table_stem)

    save_rows = True # writes to csv rater than synapse table 
    #walk_direction_preheel
    #def walk_direction_preheel(ax, ay, az, t, samplerate, 
    #                          stride_fraction=1.0/8.0, threshold=0.5,
    #                         order=4, cutoff=5, plot_test=False)
   
    

   # directions = walk_direction_preheel(ax, ay, az, t, samplerate,\
   #                                             stride_fraction, threshold0, order, cutoff)
 
    
    try:
        px, py, pz = project_walk_direction_preheel(ax, ay, az, t, samplerate,
                                                stride_fraction, threshold0, 
                                                order, cutoff)
        
        feature_row = run_pyGait(py, t, samplerate,
                                            duration, gaitthreshold, order,
                                            cutoff, distance, row, file_path,
                                            table_stem, save_rows)
    except (ValueError, IndexError) as e:
        return -1,e
    else:
        return 1,feature_row
    

                




        