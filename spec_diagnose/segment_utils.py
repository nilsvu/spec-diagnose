import os
import sys
import glob
import re
import h5py
import numpy as np


def FindLatestSegments(EvDir, Lev, tmin=-1e10, tmax=1e10, WithRingdown=True):
    """
Search through a usual SpEC segment structure, to find all
segments with start-time tmin < tastart < tmax.  Negative 
tmin count from the end of the run, e.g. use tmin=-20 to 
find the very last few segments.

if WithRingdown==False, exclude ringdown segments

Procedure:
1. Assemble list of existing directories
  [ ${EvDir}/Lev${Lev}_*] + [${EvDir}/Lev${Lev}_Ringdown/Lev${Lev}_* ]

2. Go backward through the directories
   2a. Retrieve start-time for each segment from DIR/Run/RestartTimes.txt
   2c. Terminate if start-time of current segment is
          earlier than ${DurationCompletedSegments} of start-time of last segment

3. This routine does NOT yet look for joined segments, JLev${Lev}.  Should be simple to add
JLev{$Lev} at start of list assembled in step 1., as long as CombineSegments also has combined Restart.txt


RETURNS
  segments,     -- list of *all* segments considered     [tuple of strings]
  tstart,       -- start time of each segment considered [tuple of floats]
  term_reason   -- TerminationReasons for each segment   [tuple of strings]
"""

    if not os.path.isdir(EvDir):
        raise IOError("Directory {} does not exist".format(EvDir))

    # find all segments
    tmp=os.path.join(EvDir,"Lev{}_*".format(Lev),'Run')
    all_segments=sorted(glob.glob(tmp))

    if WithRingdown:
        tmp=os.path.join(EvDir,"Lev{}_Ringdown/Lev{}_*".format(Lev,Lev),'Run')
        all_segments.extend(sorted(glob.glob(tmp)))


    segments=[]
    tstart=[]
    term_reason=[]
    tend=None
    for seg in all_segments:
        if seg[-7:]=='_AA/Run':
            # first segment (inspiral or ringdown), where RestartTimes.txt
            # is not reporting the initial start of the run
            # take Evolution.input instead
            tmp=os.path.join(seg,'Evolution.input')
            if not os.path.exists(tmp):
                raise IOError("{}--did not find Evolutiuon.input".format(seg))
            p=re.compile("^ *StartTime *= *(.+); *\n")
            for line in open(tmp):
                m=p.match(line)
                if m:
                    #print("seg={}--p.group(1)={}".format(seg,p.group(1)))
                    StartTime=float(m.group(1))
            tmp=os.path.join(seg,'RestartTimes.txt')
            if os.path.exists(tmp):
                # if we've got restarts, use them for an approximate tend
                restarts=np.loadtxt(tmp,
                                    ndmin=1   # so it always returns an array
                                    )
                tend=restarts[-1]

        else: # standard non-_AA segment
            tmp=os.path.join(seg,'RestartTimes.txt')
            if not os.path.exists(tmp):
                raise IOError("{} not found--don't yet know how to handle this".format(tmp))
            restarts=np.loadtxt(tmp,
                                ndmin=1   # so it always returns an array
                                )
            StartTime=restarts[0]
            # overwrite to find very last restart time
            tend=restarts[-1]

        tmp=os.path.join(seg,'TerminationReason.txt')
        if os.path.exists(tmp):
            with open(tmp, 'r') as myfile:
                TerminationReason = myfile.readlines()[0]
                prefix='Termination condition '
                if TerminationReason.startswith(prefix):
                    TerminationReason=TerminationReason[len(prefix):-1]
        else:
            TerminationReason='ongoing'

        segments.append(seg)
        tstart.append(StartTime)
        #print("term_reason={}".format(term_reason))
        #print("TerminationReason={}".format(TerminationReason))
        term_reason.append(TerminationReason)

    if tmin<0 and tmin!=-1e10:
        if tend is None and len(tstart)==0:
            print('specified tmin<0, which requires an estimate for tend.')
            print('No such estimate could be obtained, b/c run too short.')
            print('Therefore, tmin will be ignored.')
        else:
            tend=tstart[-1]  # use start of last segment as approx of end-time
            tmin = tend+tmin

    seg_=[]
    tstart_=[]
    term_reason_=[]
    #print("tstart={}".format(tstart))
    for s, t, r in zip(segments, tstart, term_reason):
        #print(" --- tmin={},  t={},  tmax={}".format(tmin, t, tmax))
        if t>tmin and t<tmax:
            seg_.append(s)
            tstart_.append(t)
            term_reason_.append(r)
    return seg_, tstart_,term_reason_



def LoadH5_from_segments(segments, filename, dataset_matches='',group_matches='',
                         verbose=False):
    """
Given a list of segments (incl. '/Run' directories),
check each one for a file 'filename', load that h5 file, concatenate data, and
provide it as recursive dictionary

OPTIONS:
   dataset_matches=[regex] -- only load data-sets matching the regex (e.g. 'Y_l2_m2')
                              Note:  root-datasets are always loaded
   group_matches=[regex]   -- only traverse groups matching this regex (e.g. 'R0200')
"""

    # helper function for iterative loading
    def LoadFromOpenH5(F, D,dataset_matches='',group_matches=''):
        """Given the open H5-file handle F
        (either representing the file, or a group inside the file),
        load all .dat files and store them as elements in the
        directory D.  Recursively descend into each .dir,
        and add those as sub-dictionaries into D
        dataset_matches=[regex] -- only load data-sets matching the regex"""

        for k in F.keys():
            if k.endswith('.dat') and \
            (re.match(dataset_matches,F[k].name)
             or F.parent==F # always keep top-level .dat fields
            ):
                field=k[:-4]  # remove extension from name for convenience
                data=F[k][()]
                if 'Legend' in F[k].attrs:
                    # a legend was provided for this data-set, split
                    # the dataset into individual entries, and put
                    # into a dictionary indexed by legend:
                    if not field in D:
                        D[field]={}
                    for i,legend in enumerate(F[k].attrs['Legend']):
                        legend=legend.decode("utf-8")
                        if legend in D[field]:
                            D[field][legend]=np.concatenate((D[field][legend],
                                                             data[:,[0,i]]))
                        else:
                            D[field][legend]=data[:,[0,i]]
                else:
                    # no legend for this data-set.  Load as one big array
                    if field in D:
                        D[field]=np.concatenate((D[field], data))
                    else:
                        D[field]=data
            if k.endswith('.dir') and re.match(group_matches,F[k].name):
                field=k[:-4]
                if field not in D:
                    D[field]={}
                LoadFromOpenH5(F[k],D[field],
                               dataset_matches=dataset_matches,
                               group_matches=group_matches)
        return

    D={}
    n_files=0
    for seg in segments:
        f=os.path.join(seg,filename)
        if os.path.exists(f):
            F=h5py.File(f,'r')
            LoadFromOpenH5(F,D,
                           dataset_matches=dataset_matches,
                           group_matches=group_matches)
            F.close()
            n_files=n_files+1
    if verbose:
        print("{} -- file found in {} of {} segments".format(filename,
                                                             n_files,len(segments)))
    return D



def LoadDat_from_segments_simple(segments, filename):
    """
Given a list of segments (incl. '/Run' directories),
check each one for a file 'filename', load that dat file, concatenate data, and
provide as np.array

Note: This function requires that the .dat file in all segments has
same number of columns
    """
    out=None
    n_files=0
    for seg in segments:
        f=os.path.join(seg,filename)
        if os.path.exists(f):
            tmp=np.loadtxt(f)
            if len(tmp.shape)==2: # ?? not sure why
                if out is None:
                    out=tmp
                else:
                    out=np.concatenate((out,tmp))
            n_files+=1
    print("{} -- file found in {} of {} segments".format(filename,
                                                         n_files,len(segments)))
    return out


def LoadDat_with_legend(F):
    """
Load a .dat file with standard SpEC legend strings.
Results will be placed into a dictionary indexed by the
legend string.

F -- filename

RETURNS
  D -- dictionary
"""
    if not os.path.isfile(F):
        raise IOError("File {} not found".format(F))

    # parse legend portion of file, construct legend as dictionary 'keys'
    keys={}   # dictonary of keys:  int -> legend string
    p=re.compile("^# *\[([0-9]+)\] * = *(.+)\n")
    for line in open(F):
        m=p.match(line)
        if m is not None:
            #print("{}  -{}-".format(m.group(1), m.group(2)))
            keys[ int( m.group(1) )]  = m.group(2).strip()
    #print(keys)
    tmp=np.loadtxt(F,
                   ndmin=2 # so it always returns a 2-d array
                  )
    #print(tmp.shape)
    D={}
    for idx,legend in keys.items():
        D[legend]=tmp[:,[0,idx-1]]  # -1, since SpEC legends are 1-based
    return D


def LoadDat_from_segments(segments, filename, verbose=False):
    """
Given a list of segments (incl. '/Run' directories),
check each one for a file 'filename', load that dat file, concatenate data, and
provide as np.array

The data is returned as a dictionary with 2-column data-sets, one each for each column,
named by the legend's inside the .dat files.  This allows for .dat files with a changing
number of columns, like GhCe_Linf.dat

RETURNS
  D -- dictionary
"""
    D={}

    n_files=0
    for seg in segments:
        f=os.path.join(seg,filename)
        if os.path.exists(f):
            tmp=LoadDat_with_legend(f)
            for legend, data in tmp.items():
                if legend in D:
                    D[legend]=np.concatenate((D[legend], data))
                else:
                    D[legend]=data
            n_files+=1
    if verbose:
        print("{} -- file found in {} of {} segments".format(filename,
                                                             n_files,len(segments)))
    return D


def ImportRun(path_to_ev, Lev, tmin=-1e10, tmax=1e10,verbose=False):
    """ImportRun

Load some important files for a certain Ev/Lev*, and populate a
dictionary with the imported data as follows
"""
    D={}
    segs,tstart,termination=FindLatestSegments(path_to_ev,Lev, tmin=tmin, tmax=tmax)
    D['segs']=segs
    D['tstart']=tstart
    if verbose or True:
        first_seg=segs[0]
        first_seg=first_seg[first_seg.find('Lev'):]
        last_seg=segs[-1]
        last_seg=last_seg[last_seg.find('Lev'):]
        print(f"Loading {len(segs)} segments {first_seg} @ {tstart[0]:7.3f} ... {last_seg} @ {tstart[-1]:7.3f}")
        #print("tstart={}".format(tstart))
    D['termination']=termination

    if verbose: print("Horizons",end='')
    D['Horizons']=  LoadH5_from_segments(segs,"ApparentHorizons/Horizons.h5")
    D['AhA']=       LoadDat_from_segments(segs,"ApparentHorizons/AhA.dat")
    D['AhB']=       LoadDat_from_segments(segs,"ApparentHorizons/AhB.dat")
    D['sep']=       LoadDat_from_segments(segs,"ApparentHorizons/HorizonSepMeasures.dat")
    D['ForContinuation'] = LoadDat_from_segments(segs,"ForContinuation/AhC.dat")
    if verbose: print(", GridExtents",end='')
    D['AdjustGrid']=LoadH5_from_segments(segs,"AdjustGridExtents.h5")
    if verbose: print(", Constraints",end='')
    D['GhCeLinf'] = LoadDat_from_segments(segs,"ConstraintNorms/GhCe_Linf.dat")
    if verbose: print(", DiagAhSpeeds",end='')
    D['DiagAhSpeedA'] = LoadDat_from_segments(segs,"DiagAhSpeedA.dat")
    D['DiagAhSpeedB'] = LoadDat_from_segments(segs,"DiagAhSpeedB.dat")
    if verbose: print(", DampingTimes",end='')
    D['GrAdjustMaxTstepToDampingTimes'] = \
        LoadDat_from_segments(segs,"GrAdjustMaxTstepToDampingTimes.dat")
    D['GrAdjustSubChunksToDampingTimes'] = \
        LoadDat_from_segments(segs,"GrAdjustSubChunksToDampingTimes.dat")
    D['TStepperDiag'] = LoadDat_from_segments(segs,"TStepperDiag.dat")
    return D
