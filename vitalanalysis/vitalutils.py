# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 13:50:53 2017

@author: telferm
"""
import imp
#imp.reload(fh)
#imp.reload(vs)
import pandas as pd
import numpy as np
import filehandler as fh
from scipy.signal import medfilt
def getADLwalk(walkdf,walkno,basedir,subdir,subject,samplerate=50):
    #walkno = 47
    walkdet = walkdf.iloc[walkno]
    # get filename into memory 
    fn = walkdet.fn
    sensdat = fh.getMatlabfile(fn,basedir,subdir,subject)
     #get onbody periods
             
    onbodywindow = 2 # minutes window size 
   #TODO standardise all times in seconds 
    onbodyperiods, offBodyPeriods = vs.dissectData(sensdat,onbodywindow,#
                                                       samplerate)
    #get the period of this walk (startsec, endsec within period )
    # TODO annoyingly  onbodyperiods start 1, but in walkdf they start at 0 
    # TODO even worse , the period is defined in perms of the sma window of 2 minutes 
    thisperiod = onbodyperiods[int(walkdet.period)+1]
    thisperiodstart = thisperiod[0]*60*2 + walkdet.startsec
    duration = walkdet.endsec - walkdet.startsec
    samplestart = int(thisperiodstart*samplerate)
    sampleend = samplestart + int((duration*samplerate))
     
    walkdata = sensdat[samplestart:sampleend]
    return walkdata

def getVectorsBase(subject,date,rdf2,basedir,subdir,stime='00:00:00',etime='23:59:59'):
    
    returndf = pd.DataFrame([])
    
    # id files for subject containing date requested from rdf2 
    subject = 'vah002'
    date = '2017-10-01'
    stime='00:00:00';etime='23:59:59'
    fns = rdf2[(rdf2.subject==subject) & (rdf2.date==date)]
    for i in np.arange(len(fns)):
        fn = fns.iloc[i].fn
        print('processing file:',fn)
        tempdf = fh.getMatlabfile(fn,basedir,subdir,subject) 
        returndf = pd.concat([returndf,tempdf],axis=0)
        t = tempdf.datetime.values
        t2 = np.reshape(t,(1,len(t)))
        
    #TODO file may contain more than 1 day 
    returndf.index = returndf['datetime']    
   
    returndf = returndf.truncate(before=date +' '+ stime, after = date+' '+etime)
    
    return returndf 
returndf[['x','y','z']].plot() 

    date = '2017-10-01'
    stime='11:00:00';etime='11:59:59' 
vah0021oct = vu.getVectorsBase(subject,date,rdf2,basedir,subdir,stime='00:00:00',etime='23:59:59') 
thsd,thrange =vs.getNoise(vah0021oct,subject)
vs.
#
#    returndf[['x','y','z','datetime']].plot(x='datetime')
#    tempdf3=tempdf.iloc[:20,:]  
#    plt.plot(tempdf.x)                                                  
#    # for each file: 
#    # get the data within this date (it may contain data from the day before or 
#    # the day after )  
#    #add to the df to return 
#    #onbp,offbp = onbodyperiods.copy(), offBodyPeriods.copy()
#    #imp.reload(vs)
#    onbodywindow=2
#    onbodyperiods, offBodyPeriods =  \
#        vs.dissectData(returndf,onbodywindow,#
#                       samplerate=50)
#    tempdf[['x','y','z','datetime']].plot(x='datetime')
#    tempdf2 = tempdf.copy()
#    tempdf2.index = tempdf2.datetime    
#    # looks like 7.55 - 8.00 is off body, for v001 5/10 - what is sd ? 
#    minislice = tempdf2['2017-10-05 07:55':'2017-10-05 08:00']
#    minislice[['x','y','z','datetime']].plot()  
#    minislice[['x','y','z']].std()  
#    minislice[['x','y','z']].min() 
#    minislice[['x','y','z']].max() -   minislice[['x','y','z']].min() 
#    ar = minislice[['x','y','z']].values 
#    smoothedx = medfilt(ar[:,0],3)   
#    smoothedy = medfilt(ar[:,1],3)   
#    smoothedz = medfilt(ar[:,2],3) 
#    
#       
#    #tt['sx','sy','sz'] =smoothed[:,0:2] 
#    minislice['sx'] =smoothedx
#    minislice['sy'] =smoothedy 
#    minislice['sz'] =smoothedz 
#    minislice[['sx','sy','sz']].std()
#    minislice[['sx','sy','sz']].max() - minislice[['sx','sy','sz']].min() 

def getvectors(subject,eventtype,eventnum,eventobj,data,tvalues,samplerate,offset=30):
    ''' get the acceleration vectors and times  of type:
        'raw' : not adjusted for anything
        'body' : with gravity removed 
        for a subject
        at a time - wither a walk no, pt no. or epoch 
        implementing walk no. now '''
    if eventtype=='walks':
        epoch = eventobj[eventnum]
        ssample = int((epoch[0]*samplerate)-offset) # times in walks are seconds from start 
        esample = int((epoch[1]*samplerate)+offset)
        print(esample,ssample)
        vectors = data[ssample:esample]
        tvector = tvalues[ssample:esample]
    else:
        print('I dont do that:',eventtype)
        return
            
        
    return vectors, tvector   

#%% test getADLwalk#w
#walkno = 47
#filepath = basedir+subdir
#subject = 'vah006'
#
#ondata = getADLwalk(walkdf,walkno,basedir,subdir,subject,samplerate=50)
#plt.plot(ondata[['x','y','z']])
#
#ar, arb, argr, ars, vacc, Vhacc, newdf = vs.getVectors(ondata,
#                                                       start=0,
#                                                       end=-1,
#                                                       gvalue=gravity)
#plt.plot(vacc)
#epoch = [newdf.datetime.values[0],newdf.datetime.values[-1]] #TODO these are np format - y?
#
#sma = vs.getSMA(arb)
#starttod = newdf.iloc[0].timestamp
#endtod = newdf.iloc[-1].timestamp
#t = np.linspace(starttod, endtod, len(sma))
#
#t = np.array([datetime.utcfromtimestamp(x) for x in t])
#
#actclasses=vs.getActivityClasses(sma,g=gravity) #... compute the activity classes 
#actdfs = pd.Series({'subject':subject,
#                    'fn':fn,'sma':sma,'period': periodnumber,
#                    'actclasses':actclasses,'t':t})
#actdf = actdf.append(actdfs,ignore_index=True)
#
## get walking bouts for each period of onbody time 
#
#walklog = vs.getWalks(vacc,actclasses,smawindow,epoch,samplerate=50,accelthov=1.4)
#
#
##%% test getvectors subject = 'vah006' 
#type = 'walks'
#data = vacc.copy()
##data = newdf.copy()
##data = arb.copy()
#eventnum=47
#samplerate=50
#eventobj = walks
#tvalues = newdf.timestamp.values
#vex, tv = getvectors(subject,type,eventnum,eventobj,data,tvalues,samplerate) 
#vex2 = vex[50:550]
#vex2.shape
#mm=autocorrelate(vex, unbias=2, normalize=2, plot_test=True)
#fig3 = plt.subplot(111)
#fig3.plot(vex)
#fig3.show()
#
## get horizontals 
#data = Vhacc.copy() 
#hvex, tv = getvectors(subject,type,eventnum,eventobj,data,tvalues,samplerate)   
#plt.plot(vex)  
#plt.show()     
#hmag = vs.magnitude(hvex)
#tmag = vs.magnitude(vex)
#plt.plot(hmag)
#plt.plot(tmag)
#hvex
#hvx = hvex[:,0]    
#hvy = hvex[:,1]     
#hvz = hvex[:,2] 
#ax = hvex[:,0].copy()    
#ay = hvex[:,1].copy()     
#az = hvex[:,2].copy() 
## want x to change by 1 for each step
## we have heelstrikes data at ipeaks need to incease x by 1 at each ipeak value
## so change ipeaks into a series of lengths and make x[start:end ]  of each 
## range = +1 the last 
## 3d plot of horizontal vectors for each step:
#    # 1. calculate eigenvals for each sep 
#stepper = []
#eigvecs, eigvals = [],[]
#count = 0 
#for i,x in enumerate(ipeaks):
#    toadd = [i]*(x-count)
#    cov = np.cov(hvex[count:x].T)
#    eigval,eigvec = np.linalg.eig(cov) 
#    eigvecs.append(eigvec)
#    eigvals.append(eigval)
#    stepper += toadd
#    count = x
#
#toadd = [i+1]*(len(hvx)- len(stepper))
#stepper += toadd
#    
#x = stepper
#z=0   
#y = np.arange(len(hvx))
#x = x[:126]
#y = y[:126]
#from mpl_toolkits.mplot3d import axes3d
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.quiver(x, y, z, hvx[:126], hvy[:126], hvz[:126],length = 1,  normalize=False)
#ax.quiver(0,0,0,eigx,eigy,eigz,colors=(1,0,0,1),length=1,)
#ax.set_xlabel('X axis')
#ax.set_ylabel('Y axis')
#ax.set_zlabel('Z axis')
#plt.show()
## compute eigenvalues for step ... 3 ? 
#x,y,z  = 0,0,0
#cov = np.cov(hvex.T)
#eigvals,eigvecs = np.linalg.eig(cov) 
#eigx,eigy,eigz =  eigvecs.T
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot([(0,0,0)],[(1,1,1)])
##ax.quiver(x, y, z, eigx,eigy,eigz, length =1,  normalize=False)
#ax.set_xlim(-8,8)
#ax.set_ylim(-0,126)
#ax.set_zlim(-1,1)
#ax.set_xlabel('X axis')
#ax.set_ylabel('Y axis')
#ax.set_zlabel('Z axis')
#plt.show()
