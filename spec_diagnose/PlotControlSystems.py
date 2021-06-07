import matplotlib.pyplot as plt
import numpy as np

import spec_diagnose.plot_utils as plot_utils

# from https://gist.github.com/thomasaarholt/c8440b132aaea9f71f0588af486ad457
#
# some modifications for empty axes and for log-axes
#
def autoscale(ax=None, axis='y', margin=0.1):
    '''Autoscales the x or y axis of a given matplotlib ax object
    to fit the margins set by manually limits of the other axis,
    with margins in fraction of the width of the plot
    
    Defaults to current axes object if not specified.
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    if ax is None:
        ax = plt.gca()
    newlow, newhigh = np.inf, -np.inf

    for artist in ax.collections + ax.lines:
        x,y = get_xy(artist)
        if axis == 'y':
            setlim = ax.set_ylim
            lim = ax.get_xlim()
            scale = ax.get_yscale()
            fixed, dependent = x, y
        else:
            setlim = ax.set_xlim
            lim = ax.get_ylim()
            scale = ax.get_xscale()
            fixed, dependent = y, x

        low, high = calculate_new_limit(fixed, dependent, lim)
        newlow = low if low < newlow else newlow
        newhigh = high if high > newhigh else newhigh
    
    if newlow==np.inf: 
        return  # no data at all, so don't reset anything

    # lin mapping which should be log
    if scale=='log':  
        if newlow==newhigh: # if all datapoints have identical values
            fac=1+margin
        else:
            # compute limits on log-axis to have
            # linear distance 'margin' on both ends
            fac=(newhigh/newlow)**margin
        setlim(newlow/fac, newhigh*fac)
    else:
        delta = margin*(newhigh - newlow)
        if delta==0: # if all datapoints have identical values
            delta=margin*newhigh 
        setlim(newlow-delta, newhigh+delta)
    
def calculate_new_limit(fixed, dependent, limit):
    '''Calculates the min/max of the dependent axis given 
    a fixed axis with limits
    '''
    if len(fixed) > 2:
        mask = (fixed>limit[0]) & (fixed < limit[1])
        window = dependent[mask]
        if len(window)==0:
            window=dependent # no data in plot - use all data for range
        low, high = window.min(), window.max()
    elif len(fixed)>=1:
        low = dependent[0]
        high = dependent[-1]
        if low == 0.0 and high == 1.0:
            # This is a axhline in the autoscale direction
            low = np.inf
            high = -np.inf
    else:
        # no data at all, just make up some bounds
        low=0.
        high=1.
    return low, high

def get_xy(artist):
    '''Gets the xy coordinates of a given artist
    '''
    if "Collection" in str(artist):
        x, y = artist.get_offsets().T
    elif "Line" in str(artist):
        x, y = artist.get_xdata(), artist.get_ydata()
    else:
        raise ValueError("This type of object isn't implemented yet")
    return x, y


# if __name__ == "__main__":
#     # To test
#     fig, axes = plt.subplots(ncols = 4, figsize=(12,3))
#     (ax1, ax2, ax3, ax4) = axes

#     x = np.linspace(0,100,300)
#     noise = np.random.normal(scale=0.1, size=x.shape)
#     y = 2*x + 3 + noise

#     for ax in axes:
#         ax.plot(x, y)
#         ax.scatter(x,y, color='red')
#         ax.axhline(50., ls='--', color='green')
#         for ax in axes[1:]:
#             ax.set_xlim(20,21)
#             ax.set_ylim(40,45)

#         autoscale(ax3, 'y', margin=0.1)
#         autoscale(ax4, 'x', margin=0.1)

#         ax1.set_title('Raw data')
#         ax2.set_title('Specificed limits')
#         ax3.set_title('Autoscale y')
#         ax4.set_title('Autoscale x')
#         plt.tight_layout()
#         fig.savefig('autoscale.png', dpi=400)


################################################################        
def PlotControlSystems(D, AH, tref=0., xlim=None):
    """
    Plots various diagnostics to diagnose IncomingCharFields errors.
    D -- dictionary returned form segment_utils.ImportRun()
    AH -- one of 'A', 'B'   The AH to consider
    tref -- use t-tref as x-axis.  Useful for very zoomed in time-regions late in a run
    xlim -- if given is used in call to set_xlim(xlim) to set time-axis

    
    ActivationState documented at https://github.com/sxs-collaboration/spec/blob/6e02cd7347cd34b013ebf2f0fcfe9fe85bf9acf5/Evolution/FoshSystem/DualFrameSystem/MeasureControlAhSpeed.hpp#L29
    Many other quantitites are explained at https://arxiv.org/abs/1412.1803
    """
    
    if tref==0:
        xaxis_label="t"
    else:
        xaxis_label=f"$t - {tref}$"
    

    # TOP PANELS - DIAG-AH
    fig0 = plt.figure(figsize=[32,3.5])
    diagAH=D['DiagAhSpeed'+AH]
    Panels0=[['TargetCharSpeed','Qa'],
            ['CharSpeed','ComovingMinCharSpeed','TargetCharSpeed','Qa'],
            ['CharSpeed','Q'],
            ['CharSpeedMax','ComovingMinCharSpeedMax'],
            ['DeltaR0','MinDeltaR0','ActivationState'],
            ['RelDeltaR0','MinRelDeltaR0'],
            ['DtLambdaAH','DtLambda'],
            ['LambdaAH','Lambda'],
        ]
    columns=len(Panels0)
    axs0=[fig0.add_subplot(1, columns, i+1) for i in range(columns) ]
    for k in diagAH.keys():
        found=False
        if k=='time': continue
        d=diagAH[k]
        t=d[:,0]-tref
        for idx,panel in enumerate(Panels0):
            if k in panel:
                if k in ['Lambda', 'DtLambda']:
                    axs0[idx].plot(t,-d[:,1], label="--"+k)
                else:
                    axs0[idx].plot(t,d[:,1], label=k)
                found=True
                if idx==0 and k=='Q':
                    axs0[idx].plot(t,d[:,1],'o')

        if not found: # plot in first panel
            axs0[0].plot(t,d[:,1], label=k)
    axs0[0].set_title('DiagAhSpeed'+AH)
    
    
    # MIDDLE PANELS - AH
    fig1 = plt.figure(figsize=[32,3.5]) 
    ignore1=['Area']
    Panels1=[[],
             ['NumIterations', 'convg reason'],
             ['Shape_TruncationError'],
             ['L_surface','L_mesh','L_max'],
             ['max(|r^i-c^i|)','min(|r^i-c^i|)','max(r)','min(r)'],
             ['sqrt(Area/16pi)'],
             ['Center_x'], 
             ['Center_y', 'Center_z'],
             
    ]
    columns=len(Panels1)
    axs1=[fig1.add_subplot(1, columns, i+1) for i in range(columns) ]
    for k in D['Ah'+AH].keys():
        found=False
        if k=='time' or k in ignore1: continue
        d=D['Ah'+AH][k]
        t=d[:,0]-tref
        for idx,panel in enumerate(Panels1):
            if k in panel:
                if k=='NumIterations':
                    axs1[idx].plot(t,d[:,1], '.', label=k)
                else:
                    axs1[idx].plot(t,d[:,1], label=k)

                found=True
        if not found: # plot in first panel
            axs1[0].plot(t,d[:,1], label=k)
    axs1[0].set_title('Ah'+AH)

    
    
    # BOTTOM PANELS - OTHER THINGS
    
    fig2 = plt.figure(figsize=[32,3.5]) 
    columns=8
    axs2=[fig2.add_subplot(1, columns, i+1) for i in range(columns) ]
    axs2[0].set_title('Diverse diagnostics')

    
    Panels2=[
            ['TargetChunkSize', 'ActualChunkSize'], #, 'SubChunksPerChunk'],
        ['Tdamp(LambdaFactorA)', 'Tdamp(LambdaFactorA0)', 
         'Tdamp(LambdaFactorB)', 'Tdamp(LambdaFactorB0)',],
        ['Tdamp(CutX)', 'Tdamp(ExpansionFactor)', 'Tdamp(QuatRotMatrix)', 
         'Tdamp(SkewAngle)', 'Tdamp(Trans)'],
        ['Tdamp(SmoothCoordSep)', 'Tdamp(SmoothMinDeltaRNoLam00AhA)', 'Tdamp(SmoothMinDeltaRNoLam00AhB)', 'Tdamp(SmoothRAhA)', 'Tdamp(SmoothRAhB)'],
    ]
    
    for idx,panel in enumerate(Panels2):
        for k in panel:
            d=D['GrAdjustSubChunksToDampingTimes'][k]
            t=d[:,0]-tref
            axs2[idx].plot(t,d[:,1],label=k)
    axs2[1].set_yscale('log')
    axs2[2].set_yscale('log')
    axs2[3].set_yscale('log')


    # AdjustGridExtents
    k='Sphere'+AH+'0'
    axs2[4].set_title(k)
    if k in D['AdjustGrid']:
        # === Nr ===
        # get color to use for multiple lines
        c=next(axs2[4]._get_lines.prop_cycler)['color']
        d=D['AdjustGrid'][k]['Extents']['Extent[0]']
        t=d[:,0]-tref
        axs2[4].plot(t,d[:,1]/10, color=c,
                     label="Nr/10")
        
        d=D['AdjustGrid'][k]['Bf0I1']['TruncationErrorExcess']
        t=d[:,0]-tref
        axs2[4].plot(t,d[:,1], '--', color=c,
                     label="TruncErrExcess-r")
        idx=d[:,1]>0
        axs2[4].plot(t[idx],d[idx,1], 'o', color=c)
        
        # === Ntheta ===
        # get color to use for multiple lines
        c=next(axs2[4]._get_lines.prop_cycler)['color']
        d=D['AdjustGrid'][k]['Extents']['Extent[1]']
        t=d[:,0]-tref
        axs2[4].plot(t,d[:,1]/10, color=c,
                     label="Ntheta/10")
        
        d=D['AdjustGrid'][k]['Bf1S2']['TruncationErrorExcess']
        t=d[:,0]-tref
        axs2[4].plot(t,d[:,1], '--',
                     color=c, label="TruncErrExcess-theta")
        idx=d[:,1]>0
        axs2[4].plot(t[idx],d[idx,1], 'o', color=c)

        
    # proper sep horizon
    for k in D['sep'].keys():
        if k=='t': continue
        d=D['sep'][k]
        t=d[:,0]-tref
        axs2[5].plot(t,d[:,1],label=k)

        
    # time-step size
    for k in ['dt',]:
        d=D['TStepperDiag'][k]
        t=d[:,0]-tref
        axs2[6].plot(t,d[:,1],label=k)

    for k in 'MaxAllowedTstep', 'Tdamp(LambdaFactorA0)', 'Tdamp(LambdaFactorB0)':
        d=D['GrAdjustMaxTstepToDampingTimes'][k]
        t=d[:,0]-tref
        axs2[6].plot(t,d[:,1],label=k)
    axs2[6].set_yscale('log')

    
    # common horizon in ForContinuation
    HaveContinuation = ( len(D['ForContinuation'].keys())>0 )
    axs2[7].set_title("ForContinuation/AhC")
    if HaveContinuation:
        for k in 'L_surface', 'L_mesh', 'L_max':
            d=D['ForContinuation'][k]
            t=d[:,0]-tref
            axs2[7].plot(t,d[:,1],label=k)
        tmp=axs2[7].get_xlim()
        plot_utils.AnnotateSegments(axs2[7], D, y=30, TerminationReason=True, tref=tref)
        axs2[7].set_xlim(tmp)
        axs2[7].legend(fontsize=8)  
        axs2[7].set_xlabel(xaxis_label)
        
    # OVERALL COSMETICS
    #  ForContinuation was already rescaled separately, so exclude here
    for ax in axs0+axs1+axs2[:-1]:
        ax.legend(fontsize=8)
        ax.set_xlabel(xaxis_label)
        if xlim is not None:
            ax.set_xlim(xlim)
            autoscale(ax)


