from os import chdir, mkdir, path

def mkdir_cd(dir):
    # Make a directory (if doesn't already exist) and cd to it.
    if not path.exists(dir):
        mkdir(dir)
    chdir(dir) 

exec(compile(open('check_fss.py').read(), 'check_fss.py', 'exec'))    
mkdir_cd('fss_a')
fss_a = get_fssa_data(start='2011:001')
plot_pitches(fss_a, savefig=True)
close('all')
plot_pitches_any_kalman(fss_a, savefig=True)
close('all')
mkdir_cd('../fss_b')
fss_b = get_fssb_data()
plot_pitches(fss_b, savefig=True)
close('all')
plot_pitches_any_kalman(fss_b, savefig=True)
close('all')
chdir('..')
