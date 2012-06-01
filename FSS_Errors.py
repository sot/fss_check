close('all')

msids = ['DP_PITCH', 'DP_PITCH_FSS', 'DP_ROLL', 'DP_ROLL_FSS',
         'AOSUNPRS', 'TFSSBKT1', 'TFSSBKT2', 'TCYLAFT6', 'TPC_FSSE',
         'AOFATTMD', 'AOACASEQ', 'AOPCADMD']

x = fetch.MSIDset(msids, '2008:001:00:00:00.000','2012:001:00:00:00.000') 

# filter bad times (NSM, SSM, and eclipses < 5 min)
bad_times = ['2000:048:00:00:00.000 2000:050:00:00:00.000',
             '2000:070:23:00:00.000 2000:071:13:00:00.000',
             '2001:110:00:00:00.000 2001:113:00:00:00.000',
             '2003:200:00:00:00.000 2003:201:00:00:00.000',
             '2004:315:00:00:00.000 2004:317:00:00:00.000',
             '2008:225:00:00:00.000 2008:229:00:00:00.000',
             '2008:292:00:00:00.000 2008:296:00:00:00.000',
             '2010:149:00:00:00.000 2010:153:00:00:00.000',
             '2011:187:12:00:00.000 2011:192:04:00:00.000',
             '2011:299:00:00:00.000 2011:300:12:00:00.000',
             '2002:031:16:09:58.434 2002:031:16:20:19.959',
             '2002:219:02:12:33.194 2002:219:02:22:51.644',
             '2002:229:16:49:34.857 2002:229:16:59:39.982',
             '2003:211:01:19:00.287 2003:211:01:29:09.512',
             '2003:219:00:13:37.167 2003:219:00:24:12.017',
             '2004:018:08:34:20.529 2004:018:08:44:46.154',
             '2004:034:05:51:58.039 2004:034:06:02:18.539',
             '2004:203:00:50:41.357 2004:203:01:01:27.482',
             '2006:199:07:25:17.155 2006:199:07:36:22.755',
             '2007:077:14:18:39.114 2007:077:15:02:46.139',
             '2007:011:21:08:51.203 2007:011:21:19:05.553',
             '2007:175:10:41:12.351 2007:175:10:51:14.401',
             '2007:193:23:30:34.028 2007:193:23:40:37.103',
             '2008:003:20:10:30.456 2008:003:20:22:11.931',
             '2008:003:20:12:15.006 2008:003:20:22:17.056',
             '2008:037:22:24:30.474 2008:037:22:37:37.024',
             '2008:164:19:46:43.357 2008:164:19:57:05.907',
             '2008:208:21:06:16.138 2008:208:21:47:20.713',
	     '2008:212:04:40:51.340 2008:212:05:21:23.115',
	     '2008:214:21:10:40.491 2008:214:21:51:18.416',
	     '2008:217:11:01:33.043 2008:217:11:41:36.118',
	     '2008:220:20:33:52.544 2008:220:21:14:40.719',
	     '2008:222:18:01:48.645 2008:222:18:41:55.820',
	     '2008:227:22:11:33.723 2008:227:22:51:59.348',
	     '2008:233:05:08:56.101 2008:233:05:50:36.551',
	     '2008:235:23:12:45.827 2008:235:23:53:03.252',
	     '2008:237:17:19:47.928 2008:237:18:00:13.553',
	     '2008:238:14:53:13.753 2008:238:15:33:17.853',
	     '2008:243:15:28:05.356 2008:243:16:08:09.456',
	     '2008:243:15:29:14.031 2008:243:16:09:31.456',
	     '2008:301:22:11:08.999 2008:301:22:21:13.099',
             '2008:301:22:11:14.124 2008:301:22:21:15.149',
             '2008:360:19:00:10.475 2008:360:19:10:55.575',
   	     '2009:330:13:09:27.198 2009:330:13:20:13.323',
   	     '2009:349:02:01:09.304 2009:349:02:12:17.979',
   	     '2010:044:16:05:30.073 2010:044:16:48:10.998',
	     '2011:153:01:20:24.336 2011:153:02:04:46.736',
	     '2011:187:12:23:44.974 2011:187:12:34:04.097',
   	     '2011:264:04:55:53.942 2011:264:05:38:17.442',
   	     '2005:305:12:00:00.000 2005:305:13:00:00.000',
   	     '2007:355:14:30:00.000 2007:355:15:30:00.000',   	     
   	     '2011:182:05:00:00.000 2011:182:06:00:00.000',
   	     '2003:221:16:00:00.000 2003:221:17:00:00.000',
   	     '2002:250:15:30:00.000 2002:250:17:00:00.000',
   	     '2008:037:22:00:00.000 2008:037:23:30:00.000']

for msid in msids:
    x[msid].filter_bad_times(table=bad_times)
    x[msid].filter_bad()

x.interpolate(dt=300)	

# filter for sun presence (+/- 10 min) and NPNT/STDY/KALM (for last 20 min)
sun_10min = (hstack(([False],(x['AOSUNPRS'].vals[:-1] == 'SUN '))) &
             hstack(([False, False],(x['AOSUNPRS'].vals[:-2] == 'SUN '))) &
             hstack(((x['AOSUNPRS'].vals[1:] == 'SUN '), [False])) &
             hstack(((x['AOSUNPRS'].vals[2:] == 'SUN '), [False, False])) &
             (x['AOSUNPRS'].vals == 'SUN '))
npnt_last_20min = (hstack(([False],(x['AOPCADMD'].vals[:-1] == 'NPNT'))) &
   hstack(([False, False],(x['AOPCADMD'].vals[:-2] == 'NPNT'))) &
   hstack(([False, False, False],(x['AOPCADMD'].vals[:-3] == 'NPNT'))) &
   hstack(([False, False, False, False],(x['AOPCADMD'].vals[:-4] == 'NPNT'))) &
   (x['AOPCADMD'].vals == 'NPNT'))
stdy_last_20min = (hstack(([False],(x['AOFATTMD'].vals[:-1] == 'STDY'))) &
   hstack(([False, False],(x['AOFATTMD'].vals[:-2] == 'STDY'))) &
   hstack(([False, False, False],(x['AOFATTMD'].vals[:-3] == 'STDY'))) &
   hstack(([False, False, False, False],(x['AOFATTMD'].vals[:-4] == 'STDY'))) &
   (x['AOFATTMD'].vals == 'STDY'))
kalm_last_20min = (hstack(([False],(x['AOACASEQ'].vals[:-1] == 'KALM'))) &
   hstack(([False, False],(x['AOACASEQ'].vals[:-2] == 'KALM'))) &
   hstack(([False, False, False],(x['AOACASEQ'].vals[:-3] == 'KALM'))) &
   hstack(([False, False, False, False],(x['AOACASEQ'].vals[:-4] == 'KALM'))) &
   (x['AOACASEQ'].vals == 'KALM'))

i = (sun_10min & stdy_last_20min & npnt_last_20min & kalm_last_20min)   
for msid in msids:
    x[msid].filter_bad(~i)

# define errors
pitch_error = x['DP_PITCH_FSS'].vals - x['DP_PITCH'].vals 
roll_error = x['DP_ROLL_FSS'].vals - x['DP_ROLL'].vals 

# define running standard deviation function
def running_stddev(vals,times,dt):
    t = times[0]
    i = 0
    stddev = zeros(ceil((times[-1] - times[0])/dt))
    while t < times[-1]:
        vals_i = (times > t) & (times < t + dt)
        stddev[i] = std(vals[vals_i])
        t = t + dt
        i = i + 1
    return stddev

# define daily standard deviations
pitch_error_std = running_stddev(pitch_error, x['DP_PITCH'].times, 3600*24)
roll_error_std = running_stddev(roll_error, x['DP_ROLL'].times, 3600*24)
tfssbkt1_std = running_stddev(x['TFSSBKT1'].vals, x['TFSSBKT1'].times, 3600*24)
tfssbkt2_std = running_stddev(x['TFSSBKT2'].vals, x['TFSSBKT2'].times, 3600*24)
tcylaft6_std = running_stddev(x['TCYLAFT6'].vals, x['TCYLAFT6'].times, 3600*24)
tpc_fsse_std = running_stddev(x['TPC_FSSE'].vals, x['TPC_FSSE'].times, 3600*24)

# define running average function
def running_mean(vals,times,dt):
    t = times[0]
    i = 0
    running_avg = zeros(ceil((times[-1] - times[0])/dt))
    while t < times[-1]:
        vals_i = (times > t) & (times < t + dt)
        running_avg[i] = mean(vals[vals_i])
        t = t + dt
        i = i + 1
    return running_avg

# define daily mean
pitch_error_avg = running_mean(pitch_error, x['DP_PITCH'].times, 3600*24)
roll_error_avg = running_mean(roll_error, x['DP_ROLL'].times, 3600*24)
tfssbkt1_avg = running_mean(x['TFSSBKT1'].vals, x['TFSSBKT1'].times, 3600*24)
tfssbkt2_avg = running_mean(x['TFSSBKT2'].vals, x['TFSSBKT2'].times, 3600*24)
tcylaft6_avg = running_mean(x['TCYLAFT6'].vals, x['TCYLAFT6'].times, 3600*24)
tpc_fsse_avg = running_mean(x['TPC_FSSE'].vals, x['TPC_FSSE'].times, 3600*24)

pitch_error_10avg = running_mean(pitch_error, x['DP_PITCH'].times, 3600*24*10)
roll_error_10avg = running_mean(roll_error, x['DP_ROLL'].times, 3600*24*10)

execfile('FSS_Errors_Plots.py')


#
##sun = x['AOSUNPRS'].vals
##
##entered = (sun[1:] == 'SUN ') & (sun[:-1] == 'NSUN')
##left = (sun[1:] == 'NSUN') & (sun[:-1] == 'SUN ')
##
##pitch_when_entered = x['DP_PITCH'].vals[1:][entered]
##roll_when_entered = x['DP_ROLL'].vals[1:][entered]
##t_entered = x.times[1:][entered]
##pitch_when_left = x['DP_PITCH'].vals[:-1][left]
#roll_when_left = x['DP_ROLL'].vals[:-1][left]
##t_left = x.times[:-1][left]
##
##figure(2)
##subplot(2,1,2)
##plot_cxctime(t_entered,roll_when_entered,'r.',label='roll when entered FSS FOV')
##plot_cxctime(t_left,roll_when_left,'b.',label='roll when left FSS FOV')
##title('Roll Angle At FSS FOV Boundary')
##ylabel('Roll [deg]')
##subplot(2,1,1)
##plot_cxctime(t_entered,pitch_when_entered,'r.',label='pitch when entered FSS FOV')
##plot_cxctime(t_left,pitch_when_left,'b.',label='pitch when left FSS FOV')
##title('Pitch Angle At FSS FOV Boundary')
##ylabel('Pitch [deg]')
##
##figure(3)
##plot(roll_when_entered,pitch_when_entered,'r.')
##plot(roll_when_left,pitch_when_left,'b.')
##title('Pitch and Roll Angle at FSS FOV Boundary')
##xlabel('Roll [deg]')
##ylabel('Pitch [deg]')
#
