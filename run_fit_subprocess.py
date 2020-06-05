'''
The reflection fit done automaticallyd
This will only work with the central run
'''

import os, sys, pickle, iminuit, time, subprocess
import numpy as np
from scipy import optimize
from copy import deepcopy
sys.path.append('/home/jpyanez/snoplus/snoplus_python')
import soc_tools

linear_fit = True
python_minimizer = True


if len(sys.argv) > 1:
    debug = True
else:
    debug = False

# Get a some unique number for the files and logs
if debug:
    unique_id = '8811'
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
#input_macro      = '/home/jpyanez/snoplus/analysis/greydisc_reflections/reflection_fit/laserball_115225_template.mac' ## 300 sec for 200 flashes
#input_macro      = '/home/jpyanez/snoplus/analysis/greydisc_reflections/reflection_fit/laserball_115225_template_simple_v1.mac' # 238 sec for 200 flashes
input_macro      = '/home/jpyanez/snoplus/analysis/greydisc_reflections/reflection_fit/laserball_115225_template_simple_v2.mac' # 252 for 200 flashes
output_macro     = os.path.join(base_dir, 'temp_macro.mac')
temp_socfile     = os.path.join(base_dir, 'temp_socfile.soc')
temp_pickle_file = os.path.join(base_dir, 'temp_toa.pckl')
temp_results     = os.path.join(base_dir, 'results.pckl')
mylog            = os.path.join(base_dir, 'mylog.log')
toa_script       = '/home/jpyanez/snoplus/analysis/greydisc_reflections/reflection_fit/get_toa.py'

nflashes         = '50'
run_nr= '115225'



# Try again
# Parameters of the fit
parameters = {
    #'spec_proba_m': 0.23,
    'spec_proba_b': 0.20, #0.13,
    'spec_std':    0.10, #, 0.24,
    'nonspec_std': 0.40, #0.34, #0.4,
    'nonspec_std_focused': 0.001, #0.00976,
    'nonspec_focused_prob': 0.30, #0.15,
    'nonspec_b':  0.5, #0.67
    }
params = parameters.keys()

kwargs = kwargs = {'forced_parameters': params,
                   'use_array_call':True,
                   'errordef':1.}

pvalues = []
bounds  = []
for p in params:
    kwargs[p] = parameters[p]
    kwargs['error_'+p] = 1.0
    kwargs['limit_'+p] = (0.0001, 0.8)
    pvalues.append(parameters[p])
    bounds.append([0.001, 0.8])
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
chi2bool = (data['time_edges'][:-1]>15.)*(data['time_edges'][:-1]<80.)
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
                this_string = base_string + one_key + ' ' + "%0.7f" % params[one_key] + '\n'
                outfile.write(this_string)   
        else:
            if default_gd_table in line:
                line = line.replace(default_gd_table, gd_table)
                #print line
            outfile.write(line)
    outfile.close()
    infile.close()


def run_macro():
    macro_string = ['rat',
                    output_macro,
                    '-P','-X',
                    '-n', run_nr,
                    '-N', nflashes,
                    '-o', temp_socfile,
                    '-l', temp_socfile+'.log']
    print 'Running macro '
    print macro_string
    subprocess.check_call(macro_string)
    
    

def make_toa():
    if False:
        soc_tools.getTOAhistogram_detailed(infile_list = [temp_socfile],
                                           outfile=temp_pickle_file)
    toa_string = ['python',
                  toa_script,
                  temp_socfile,
                  temp_pickle_file]
    subprocess.check_call(toa_string)

def chi2():
    mc_new = soc_tools.loadTOAhistogram_detailed(temp_pickle_file)
    mc_toa = mc_new['toa'].sum(axis=0)

    approx_scale = mc_toa[chi2bool].sum()/data_toa[chi2bool].sum()

    chi2_sum = 1E9
    scale_used = 0

    # Printing data used
    #print data_toa[chi2bool]
    #print data_toa[chi2bool]*approx_scale
    #print mc_toa[chi2bool]
    
    if np.sum(mc_toa[chi2bool] == 0) >= 1:
        print 'Problem with chi2, not making enough events!!'
        print mc_tao[chi2bool]
        return 1E9
    
    for test_scale in np.linspace(0.8, 1.2, 1000):

        # Linear
        if linear_fit:
            difference = data_toa[chi2bool]*approx_scale*test_scale - mc_toa[chi2bool]
            error2 = mc_toa[chi2bool]

        # Logarithmic
        else:
            difference = np.log10(data_toa[chi2bool]*approx_scale*test_scale)-np.log10(mc_toa[chi2bool])
            error2 = np.log10(mc_toa[chi2bool])
        
        
        test_chi2 = difference**2/error2
        
        if test_chi2.sum() < chi2_sum:
            chi2_sum  = test_chi2.sum()
            chi2_full = deepcopy(test_chi2)
            scale_used = deepcopy(test_scale)

    print '******** SCALE USED: ', scale_used, '******'
    return chi2_sum

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
            print 'Why is this a NAN? ', one_key
            print 'Returning a really large likelihood value'
            return 1E9
        param_set[one_key] = test_values[ikey]

        
    make_macro(param_set)
    if not debug:
        run_macro()
        #time.sleep(10)
        make_toa()
        #try:
    #make_toa()
        #except:
        #    print 'TOA did not work'

    if not os.path.isfile(temp_pickle_file):
        for i in range(20):
            print '\n\n\n\n\n******The TOA file is missing!************'
        raw_input()
        make_toa()

    try:
        os.remove(temp_socfile)
    except:
        print 'No SOC file to remove'
        
    this_chi2 = chi2()
    logstring2 = '\t\t'.join(["%0.7f" % this_chi2] + ["%0.7f" % x for x in test_values])
    
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

# For testing
#fcn([0.1, 0.2, 0.3, 0.01, 0.1, 0.6])

# This is the iminuit way

if python_minimizer:
    m = iminuit.Minuit(fcn,
                       **kwargs)
    r = m.migrad(ncall=500)
    pickle.dump({'result':r,
                 'arguments':kwargs},
                open(temp_results,'w'))

else:
    r = optimize.minimize(fcn, x0 = pvalues,
                          method = 'L-BFGS-B', # TNC, SLSQP, COBYLA
                          bounds = bounds)
    pickle.dump({'result':r,
                 'arguments':kwargs},
                open(temp_results,'w'))
                          
    
logfile.close()

### Again, rewriting the code
## First make the macro
## Then run the macro
## Then produce the pickle file
## Load pickle file, compare with data
