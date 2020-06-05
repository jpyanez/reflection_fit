'''
The reflection fit done automaticallyd
This will only work with the central run
'''

import os, sys, pickle, iminuit, time, subprocess
import numpy as np
sys.path.append('/home/jpyanez/snoplus/snoplus_python')
import soc_tools


if len(sys.argv) > 1:
    debug = True
else:
    debug = False

# Get a some unique number for the files and logs
if debug:
    unique_id = '666'
else:
    unique_id = "%i" % (10000*np.random.rand())
base_dir = os.path.join('/home/jpyanez/scratch/temp',unique_id)

if not debug:
    if os.path.isdir(base_dir):
        print 'Error!'
        sys.exit()
    else:
        os.mkdir(base_dir)


# Horrible global variables
default_gd_table = 'DiscOpticsTuned_InAV'
#gd_table         = 'DiscOpticsTuned_InAV' #'TuneOptics_31'
gd_table         = 'TuneOptics_32' # This is -2% ish
input_macro      = 'laserball_115225_template.mac'
output_macro     = os.path.join(base_dir, 'temp_macro.mac')
temp_socfile     = os.path.join(base_dir, 'temp_socfile.soc')
temp_pickle_file = os.path.join(base_dir, 'temp_toa.pckl')
temp_results     = os.path.join(base_dir, 'results.pckl')
mylog            = os.path.join(base_dir, 'mylog.log')

nflashes         = '10'
run_nr= '115225'



# Try again
# Parameters of the fit
parameters = {
    #'spec_proba_m': 0.23,
    'spec_proba_b': 0.13,
    'spec_std':    0.24,
    'nonspec_std': 0.34, #0.4,
    'nonspec_std_focused': 0.01,
    'nonspec_focused_prob': 0.15,
    'nonspec_b':  0.67
    }
params = parameters.keys()
kwargs = kwargs = {'forced_parameters': params,
                   'use_array_call':True,
                   'errordef':1.}
for p in params:
    kwargs[p] = parameters[p]
    kwargs['error_'+p] = 0.05
    kwargs['limit_'+p] = (0.0001, 0.8)
    #print params


### Data information for chi2 calculation
data_dir = '/home/jpyanez/snoplus/data/snoplus_data/laserball_runs/soc_histograms_detailed'
infile_name = os.path.join(data_dir,  run_nr+ '_TOA.pckl')
data = soc_tools.loadTOAhistogram_detailed(infile_name)
data_toa = data['toa'].sum(axis=0)

# Using the late pulsing to match
dbool = (data['time_edges'][:-1]>15.)*(data['time_edges'][:-1]<35.)
#dbool = (data['time_edges'][:-1]>-5.)*(data['time_edges'][:-1]<10.)

# Chi-squared region for comparison
chi2bool = (data['time_edges'][:-1]>35.)*(data['time_edges'][:-1]<80.)
scale = data_toa[dbool].sum()


def make_macro(params):
    
    infile = open(input_macro)
    outfile = open(output_macro, 'w')
    
    for line in infile:
        if 'PARAMETERS_TO_TUNE' in line:
            outfile.write(line)
            # Write all parameters
            base_string = '/rat/db/set GREY_DISC_PARAMETERS['+gd_table+'] '
            for one_key in params:
                this_string = base_string + one_key + ' ' + "%0.5f" % params[one_key] + '\n'
                outfile.write(this_string)   
        else:
            if default_gd_table in line:
                line = line.replace(default_gd_table, gd_table)
                #print line
            outfile.write(line)


def run_macro():
    macro_string = ' '.join(['rat ',
                    output_macro,
                    ' -P -X',
                    ' -n '+ run_nr,
                    ' -N '+ nflashes,
                    ' -o '+ temp_socfile,
                   ' -l '+ temp_socfile+'.log'])
    print 'Running macro '
    print macro_string
    os.system(macro_string)
    

def make_toa():
    soc_tools.getTOAhistogram_detailed(infile_list = [temp_socfile],
                                       outfile=temp_pickle_file)
    

def chi2():
    mc_new = soc_tools.loadTOAhistogram_detailed(temp_pickle_file)
    mc_toa = mc_new['toa'].sum(axis=0)

    this_scale = mc_toa[dbool].sum()*1./scale
    chi2_full = (data_toa[chi2bool]*this_scale - mc_toa[chi2bool])**2/ mc_toa[chi2bool]

    if np.sum(np.isnan(chi2_full))>2:
        print 'Problem with chi2!!!'
        print chi2_full
        
    
    return np.sum(chi2_full[~np.isnan(chi2_full)])

#do_header = True

def fcn(test_values):
    global do_header
    param_set = {}
    print '\nNew evaluation ', test_values

    logstring1 = '\t'.join(['fval']+params)
    logaux = ' '.join(["%0.3f" % x for x in test_values])
    
    for ikey, one_key in enumerate(params):
        if np.isnan(test_values[ikey]):
            print '\n\n************* ERROR!!!! *************'
            print 'Why is this a NAN?'
            print 'Returning a really large likelihood value'
            return 1E9
        param_set[one_key] = test_values[ikey]

        
    make_macro(param_set)
    if not debug:
        run_macro()
        time.sleep(10)
        make_toa()
        ## Delete the temporary soc file
        os.remove(temp_socfile)

    if not os.path.isfile(temp_pickle_file):
        for i in range(20):
            print '\n\n\n\n\n******The TOA file is missing!************'
        raw_input()
        make_toa()
        
    this_chi2 = chi2()
    logstring2 = '\t\t'.join(["%0.5f" % this_chi2] + ["%0.5f" % x for x in test_values])
    
    logfile = open(mylog,'a')
    logfile.write(logstring2+'\n')
    logfile.close()
    print logstring1
    print logstring2
    #print 'Chi squared ', this_chi2
    return this_chi2

logfile = open(mylog, 'w')
logstring1 = '\t'.join(['fval']+params)
logfile.write(logstring1+'\n')
logfile.close()



# This is the iminuit way
m = iminuit.Minuit(fcn,
                   **kwargs)
r = m.migrad(ncall=500)
pickle.dump({'result':r,
             'arguments':kwargs},
            open(temp_results,'w'))

logfile.close()

### Again, rewriting the code
## First make the macro
## Then run the macro
## Then produce the pickle file
## Load pickle file, compare with data
