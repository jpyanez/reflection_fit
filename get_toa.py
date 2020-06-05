import os, sys, pickle
import numpy as np


sys.path.append('/home/jpyanez/snoplus/snoplus_python')

import soc_tools

socfile_name = sys.argv[1]
outfile_name = sys.argv[2]

print 'TOA extraction'
print socfile_name
print outfile_name

soc_tools.getTOAhistogram_detailed(infile_list = [socfile_name],
                                   outfile=outfile_name)
