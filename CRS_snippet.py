#Last update: May 14, 2023: adapt to CIRES job interview 
#             January 12, 2023: Use corrected velocity to automatically find legs
#             December 20, 2022: add MaskCoPol
#             October 19, 2022: add scaledown to scale down airborne data to spaceborne
#             August 17, 2022: 
def read(filename=''):
    """
    Read and calcualte the altitude for each bin of each vertical profile.
    Radar_range is the distance from aircraft to the each radar bin, first bin is closest to aircraft, about 5000 meters.
    Convert lon, lat, and timeutc to 2d for plotting, though lon, lat, timeutc donnot change vertically.
    """
    import h5py
    import numpy as np
    import pandas as pd

    with h5py.File(filename, 'r') as f:
        #list(f.keys())
        aircraft = f['Information/Aircraft']     # NASA ER-2
        radarname = f['Information/RadarName']   # 
        datacontact = f['Information/DataContact']
        date = f['Information/FlightDate']

        lat = np.array(f['Navigation/Data/Latitude'])
        lon = np.array(f['Navigation/Data/Longitude'])
        height = np.array(f['Navigation/Data/Height'])  # Aircraft height above sea level in meters, float32
        track = f['Navigation/Data/Track']   # Direction from motion in degrees from north. 90 degrees is eastward motion

        dBZ = np.array(f['Products/Data/dBZe'])
        vel_cor = np.array(f['Products/Data/Velocity_corrected'])
        vel_uncor = np.array(f['Products/Data/Velocity_uncorrected'])
        swidth = np.array(f['Products/Data/SpectrumWidth'])
        ldr = np.array(f['Products/Data/LDR'])
        sigma0 = np.array(f['Products/Data/sigma0'])

        radar_range = np.array(f['Products/Information/Range'])    # Range in meters from the aircraft, float32
        gate_spacing = f['/Products/Information/GateSpacing']  # It is a constant for the file I checked (HIWRAP, Feb01, 26.25 meters.
        gate_spacing_des = f['/Products/Information/GateSpacing_description']
        gate_spacing_unit = f['/Products/Information/GateSpacing_units']
        maskCoPol = f['/Products/Information/MaskCoPol'][()]
        maskCrPol = f['/Products/Information/MaskCrPol'][()]

        pri = f['/Products/Information/PRI']
        #pri_des = f['/Products/Information/PRI_description'] #the 2020 data final version (RevB) does not have PRI_description 

        timeutc = np.array(f['Time/Data/TimeUTC'])  # UTC profile time in Unix Epoch format (seconds since 1970). 
                                          # Obtained from aircraft NTP

    nrange = radar_range.size
    ntime = height.size
    print('ntime = ', ntime, ' nrange = ', nrange)
    bin_alt = np.empty((ntime,nrange))
    lon_2d = np.empty((ntime,nrange))
    lat_2d = np.empty((ntime,nrange))
    timeutc_2d = np.empty((ntime,nrange))
    for itime in range(ntime):
        bin_alt[itime,:] = height[itime] - radar_range[:] # bin_altitude = Aircraft height - radar range
    for irange in range (nrange):
        lon_2d[:,irange] = np.squeeze(lon[:])
        lat_2d[:,irange] = np.squeeze(lat[:])
        timeutc_2d[:,irange] = np.squeeze(timeutc[:])

    dti = pd.to_datetime(timeutc,unit='s').round('250 ms')

    return dBZ,vel_cor,vel_uncor,swidth,ldr,timeutc,dti,timeutc_2d,lat_2d,lon_2d,radar_range,height,bin_alt,maskCoPol,maskCrPol

def info(timeutc, info=None):
    """
    This function is to understand how time is recorded in the data:
    """
    
    import numpy as np
    import datetime as dt

    # Try to use datetime module to check the start time and end time of this h5 file.
    # I may just use simpler modules like calendar and time, instead of datetime
    if info:
        print('Record starts at {} Epoch time, corresponding to {}'
              .format(timeutc[0],dt.datetime.utcfromtimestamp(timeutc[0])))
        print('Record ends at {} Epoch time, corresponding to {}'
              .format(timeutc[-1],dt.datetime.utcfromtimestamp(timeutc[len(timeutc)-1])))

    # Since I don't know if the vertical profiles are evenly distributed in time, I calculate the 
    # time differences, dt, and found all the time difference between profils are 0.5 seconds.
        print()
        print('First and 2nd records are at: ', float(timeutc[0]),float(timeutc[1])) # Use float function, so that it print out the longest number of digits, 
                                                   # otherwise, it only output a truncated number of digits, which can not
                                                   # show the time interval between two profile is 0.5 seconds.
        n = len(timeutc)
        delta_t = np.zeros(n-1,dtype='float')
        for i in range(n-1):
            delta_t[i] = timeutc[i+1] - timeutc[i]

        print('The number of time intervals that are equal to 0.25 seconds: ',np.count_nonzero(delta_t == 0.25)) # print out the number of values equal to 0.25 in array dt 
        print('Therefore, we have confirmed that all the profile records are evenly recorded at 0.25 seconds interval.')
        print('Indeed, I\'ve already rounded pandas dti (DatetimeIndex) to 0.25 seconds and output from the read function.') 
    
    return

def auto_leg(vel_cor,dti,savetxt=None):

    """
    Use corrected Doppler velocity to automatically find legs. i.e. where not all data in a column is NaN
    But it is limited to how the Doppler velocity is corrected. 
    It may chunk out some legs that I would look at a longer distance. E.g., for the Feb 8, 2022 case, GPM underflight leg.
    """
    
    import pandas as pd
    from operator import attrgetter

    df = pd.DataFrame(vel_cor)
    sr = df.any(axis=1) # For each column (time step), if any non-NaN value exist, give True to that time step, and 
                       # put the assigned True or False, for each time step, to a Series, sr

    grp = sr.eq(False).cumsum() #If a group of True is connected to each other, give a same digital number to them

    idx_s_e = grp.loc[sr.eq(True)].groupby(grp).apply(lambda x: [x.index.min(),x.index.max()])  #start and end index of all the groups
    dti_s_e = []
    dti_s_e_tuple = []

    attrs = ('year','month','day','hour','minute','second','microsecond')

    print('{} legs are automatically identified based on corrected Doppler velocity.'.format(len(idx_s_e)))
    for i in idx_s_e.index:
        #print out the datetime:
        #print(idx_s_e[i], '  ',dti[idx_s_e[i][0]],dti[idx_s_e[i][1]])
        #print out the tuple from datetime:
        print(idx_s_e[i], '  ',attrgetter(*attrs)(dti[idx_s_e[i][0]]),attrgetter(*attrs)(dti[idx_s_e[i][1]]))
        dti_s_e.append([dti[idx_s_e[i][0]],dti[idx_s_e[i][1]]])
        dti_s_e_tuple.append([attrgetter(*attrs)(dti[idx_s_e[i][0]]),attrgetter(*attrs)(dti[idx_s_e[i][1]])])
    print('Need to make sure that the automatically identified legs make sense.')

    CRS_date = dti_s_e[0][0].strftime('%Y%m%d')
    num_leg = len(idx_s_e)
    if savetxt:
        with open(savetxt+CRS_date+'_CRS_auto_leg.txt','w') as f:
           f.write(CRS_date+': '+str(num_leg)+' CRS auto legs are identified based on corrected Doppler velocity: \n')
           for i in idx_s_e.index:
               #f.writelines(str(idx_s_e[i])+'\n') #, '  ',attrgetter(*attrs)(dti[idx_s_e[i][0]]),attrgetter(*attrs)(dti[idx_s_e[i][1]]))
               f.writelines(str(idx_s_e[i])+'  '+str(attrgetter(*attrs)(dti[idx_s_e[i][0]]))+' '+str(attrgetter(*attrs)(dti[idx_s_e[i][1]]))+'\n')
           f.write('Note: Make sure the above start and end times of legs make sense. It may be over restrictive.\n')
           f.write('      To find the corresponding CPL leg, note that \n')
           f.write('      CPL has a record every 5 seconds (Level-2), or 1 seconds (Level-1), just generally.')

    return list(idx_s_e), dti_s_e, dti_s_e_tuple   #return list

def time_range(start_time,end_time,timeutc=None,note=None):
    """
    Find the start and end index
    """
    # Convert the tuple start_time and end_time to Epoch utc time
    
    import calendar
    import numpy as np
    import datetime as dt
    
    s_timeutc = calendar.timegm(start_time)
    e_timeutc = calendar.timegm(end_time)

    s = int(np.where(timeutc == s_timeutc)[0])
    e = int(np.where(timeutc == e_timeutc)[0])
    start_time_str = str(dt.datetime.utcfromtimestamp(s_timeutc))
    end_time_str = str(dt.datetime.utcfromtimestamp(e_timeutc))
    if note == 'more':
        print('The start and end time in decimal seconds since the epoch are ',timeutc[s],timeutc[e])
    print('Start at {}     end at {}'.format(start_time_str,end_time_str))
    if note == 'more':
        print('The start (s) and end (e) indices are:',s,e)

    return s,e,start_time_str,end_time_str

def plot(timeutc,timeutc_2d,lat_2d,lon_2d,bin_alt,s,e,start_time_str,end_time_str,
         var=None,
         savefig=None,fig_path='',P3=None,P3_tick=None,P3_tick_special=None,sonde=None,model=None,TAMMS=None,
         CPL=None,
         xaxis='Time',timetick=None,lattick=None,lontick=None,fs=14,man_title=None,**param_dict):
    """
    plot one variable: var, it can be dBZ, vel_cor, and other variables in the CRS data
    """

    import numpy as np
    import time
    import datetime as dt
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.ticker import MaxNLocator
    from matplotlib.ticker import MultipleLocator
    from matplotlib.ticker import AutoMinorLocator
    from matplotlib import colors

    c = var[0][s:e+1,]
    c_min = np.nanmin(c)
    c_max = np.nanmax(c)

    if var[1][0] == 'fix':
        levels = np.linspace(var[1][1][0],var[1][1][1],var[1][1][2])
    if var[1][0] == 'auto':
        print('Automatically select color levels.')
        levels = MaxNLocator(nbins=30).tick_values(c_min, c_max) #this will include c_min and c_max within the levels

    cmap = plt.get_cmap('nipy_spectral')
    norm = colors.BoundaryNorm(levels, ncolors=cmap.N)

    if 'Doppler' in var[3]:
        #cmap='seismic'
        #cmap='coolwarm'
        cmap=var[1][2] # choose 'seismic' or 'coolwarm' in the notebook
        if var[1][0] == 'fix':
            min_level=var[1][1][0]
            max_level=var[1][1][1]
        if var[1][0] == 'auto':
            print('Automatically select color levels.')
            min_level=c_min
            max_level=c_max
        divnorm=colors.TwoSlopeNorm(vmin=min_level, vcenter=0., vmax=max_level)
        norm=divnorm

    if 'Mask' in var[2]:
        levels = np.linspace(-0.5,3.5,5) #give 5 levels that make the cloud phase 0,1,2,3 in the middle of each bin
        cmap = colors.ListedColormap(['white','dimgray','darkgray','lightgray'])
        norm = colors.BoundaryNorm(levels, ncolors=cmap.N)
        
    fig = plt.figure(figsize=(15,5))
    ax = plt.axes()
    if xaxis == 'Index':
        index1d = np.arange(c.shape[0]) #this index begins at the start of the leg
        index2d = np.zeros((c.shape[0],c.shape[1]))
        for i in range(c.shape[1]):
          index2d[:,i] = index1d
        panel = ax.pcolormesh(index2d,bin_alt[s:e+1,]/1000,c,cmap=cmap,norm=norm, shading='auto')
        ax.set_xlabel('Index',fontsize=fs)
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    if xaxis == 'Time':
        if time.strftime('%S',time.gmtime(timeutc[s])) == '00':
            print('Start at a whole minute.')
            min_tick = np.arange(timeutc[s],timeutc[e],60)
        else:
            ......
        min_tick_label = []
        for i in min_tick:
            min_tick_label.append(time.strftime('%H:%M',time.gmtime(i)))
        panel = ax.pcolormesh(timeutc_2d[s:e+1,],bin_alt[s:e+1,]/1000,c,cmap=cmap,norm=norm, shading='auto')
        ax.set_xlabel('Time (UTC)',fontsize=fs)
        if timetick:
            ax.set_xticks(min_tick[::timetick[0]]) # frequency to give time tick, default 4 min
            ax.set_xticklabels(min_tick_label[::timetick[0]]) # frequency to label time tick, default 4 min
            ax.xaxis.set_minor_locator(MultipleLocator(timetick[1])) # set minor ticks to be every timetick[1], e.g.,10 seconds
        else:
            ax.set_xticks(min_tick[::5])
            ax.set_xticklabels(min_tick_label[::5])
            ax.xaxis.set_minor_locator(AutoMinorLocator())
    if xaxis == 'Lat':
        .......
    if xaxis == 'Lon':
        .......

    ax.tick_params(axis='both', which='major', labelsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)

    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylim(bottom=0,top=12) 
    ax.set_ylabel('Height (km)',fontsize=fs)

    ax.set_title('CRS '+ var[2] +' ' + start_time_str[0:10] \
              + ' ' + start_time_str[11:19] + ' to ' + end_time_str[11:19] + ' UTC' + '  ' + var[6])
    if man_title:
       ax.set_title(man_title[0],fontsize=man_title[1])
    
    if var[5][0] == 'maxmin': # note the max and min when pass the string 'maxmin'
        if var[5][1] == 'in':
            ax.text(0.03,0.93,'max = '+str(round(c_max,2)),transform=ax.transAxes)
            ax.text(0.03,0.88,'min = '+str(round(c_min,2)),transform=ax.transAxes)
        if var[5][1] == 'out':
            ax.text(0.05,-0.108,'min = '+str(round(c_min,2)),transform=ax.transAxes)
            ax.text(0.85,-0.108,'max = '+str(round(c_max,2)),transform=ax.transAxes)

    if TAMMS: #only works when xaxis='Lat'
        ......

    if CPL:
        ......

    if P3: #only can overlay P3 track when xaxis == 'Lat' or 'Lon'
        ......
        if P3['note']: #Add hh:mm:ss on the vertical cross section
            if xaxis == 'Lat':
                ......
            if xaxis == 'Lon':
                ......
            if xaxis == 'Time':
                ......
    if sonde and xaxis == 'Lat':
        ......
    if P3_tick:
        ......

    if model:
        cross = model['HRRR']
        lat2d = np.vstack([cross.latitude]*40)
        if model['allT']: 
            cl_2 = ax.contour(lat2d,cross.gh/1000., cross.t-273.15, colors='k',
                 levels=np.concatenate((np.arange(-60,0,5),np.arange(5,20,5))),
                      linestyles='dashed',linewidths=1)
            ax.clabel(cl_2, fmt='%1.0f')
        cl_3 = ax.contour(lat2d,cross.gh/1000., cross.t-273.15, colors='r',
                  levels=[0],linestyles='dashed',linewidths=2)
        ax.clabel(cl_3, fmt='%1.0f')

    ax.grid(True,**param_dict)
    #fig_name = fig_name_case + '_' + start_time_str[11:13] + start_time_str[14:16] \
    fig_name = 'CRS_'+var[4]+'_'+start_time_str[11:13]+start_time_str[14:16]+start_time_str[17:20] \
               +'-'+end_time_str[11:13]+end_time_str[14:16]+end_time_str[17:20]+ 'UTC'+'_'+xaxis

    if savefig:
        plt.savefig(fig_path+fig_name+savefig, dpi=300,bbox_inches='tight')

def scaledown(s,e,hstep,vstep,bin_alt,lat_2d,var=None,sens=None):
    ......


