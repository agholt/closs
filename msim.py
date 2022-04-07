#!/usr/bin/env python3

__author__ = "Alan Holt"
__copyright__ = "Copyright 2019"
__credits__ = ["Adrian Davies"]
__license__ = "Proprietary"
__version__ = "1.4.0"
__maintainer__ = "Alan Holt"
__email__ = "agholt@gmail.com"
__status__ = "Production"


#from Bio.Seq import Seq
#from Bio import SeqIO
#from Bio import SeqUtils
#from Bio import Entrez
#from Bio.Alphabet import IUPAC
#from Bio.Data import CodonTable


import itertools
import datetime
import json
from io import StringIO
import hashlib
import timeit
import csv
import sys


from optparse import OptionParser


import random
import numpy as np

import Msim as msim


def proc_options():

    parser = OptionParser()
    parser.add_option('-c', '--count', dest='count',
                      help = 'specify number of repetitions')
    parser.add_option('-d', '--desc', dest='desc',
                      help = 'specify description')
    parser.add_option('-p', '--params', dest='paramfile',
                      help = 'specify parameter file')
    parser.add_option('-r', '--runlength', dest='runlength',
                      help = 'runlength')
    parser.add_option('-D', '--Debug', dest='debug',
                      help = 'debug level')
    parser.add_option('-m', '--message', dest='message',
                      help = 'message')

    return parser.parse_args(sys.argv[1:])




if __name__== '__main__':


    #print("Author: %s"  % (__author__))
    #print("Credits: %s" % (__credits__))
    #print("Version: %s" % (__version__))
    #print("Status: %s"   % (__status__))
    #print("Copyright: %s" % (__copyright__))


    options,files = proc_options()

    if options.paramfile != None:
        msim_file = options.paramfile
    else:
        msim_file = "msim.json"

    if options.count != None:
        count = int(options.count)
    else:
        count = 1

    if options.desc != None:
        desc = options.desc
    else:
        desc = "test"


    if options.runlength != None:
        runlength = options.runlength
    else:
        runlength = None


    if options.debug != None:
        msim.DEBUG = int(options.debug)
    else:
        msim.DEBUG = int(0)

    if options.message != None:
        print(options.message) 


    msg = "DEBUG: %d" % (msim.DEBUG)
    msim.write_debug(msg, 0)

    msg = "command-line args: %s %s %s" % (msim_file, count, desc)
    msim.write_debug(msg, 1)



    # Get paramaters from parameter file
    with open(msim_file, 'rU') as fd:
        params = json.load(fd)
    meta_data = params['META']
    mito_params = params['MITO_PARAMS']
    #mtdna_params = params['MTDNA_PARAMS']

    msg = "%s" % (meta_data)
    msim.write_debug(msg, 1)
    msg = "%s" % (mito_params)
    msim.write_debug(msg, 1)

    results_dir = "."

    sim_time_params = meta_data['sim_time_params']

    sample_per_hour = int(sim_time_params['sample_per_hour'])
    no_days = int(sim_time_params['no_days'])
    no_years = int(sim_time_params['no_years'])


    # if simulation run length specied as option, override run length from param file
    if not runlength:
        run_length = int(sim_time_params['run_length'])
    else: 
        run_length = int(runlength)

    msg = "simulation run legnth is %s" % (run_length)
    msim.write_debug(msg, 2)

    sim_time_params = meta_data['sim_time_params']
    sample_per_hour = int(sim_time_params['sample_per_hour'])
    no_days = int(sim_time_params['no_days'])
    no_years = int(sim_time_params['no_years'])
    run_length = int(sim_time_params['run_length'])

    wild_seq = mito_params['MTDNA_PARAMS']['WILD_SEQ']
    wild_len = len(wild_seq); wild_len
    wild_type = hashlib.sha224(wild_seq.encode('utf-8')).hexdigest()


    init_pop_size = meta_data['INIT_POP_SIZE']

    msg = "wild_type: %s (%s)" % (wild_type, wild_len)
    msim.write_debug(msg, 1)

    msg = "Results directory: %s" % (results_dir)
    msim.write_debug(msg, 2)


    # Output files
    ts = str(datetime.datetime.now()).split('.')[0].replace(' ', '_')
    log_file = "msim-%s_%s.log" % (desc, ts); log_file
    mito_params['LOGFILE'] = log_file


    results_file = "%s/results-%s_%s.json" % (results_dir,desc,ts)

    # Simulation ID    
    sim_id = "%s-%s" % (desc,ts)

    results_lst = []

    for i in range(count):

        msg = "Simulation run: %d" % (i)
        msim.write_debug(msg, 1)

        # Initialise
        ts_file = "%s/ts-%s_%s_%s.json" % (results_dir,desc,ts,i)
        gen_file = "%s/gen-%s_%s_%s.json" % (results_dir,desc,ts,i)
        state_file = "%s/state-%s_%s_%s.json" % (results_dir,desc,ts,i)


        O = msim.Mito(meta_data,mito_params)
        #seq_lst = [(wild_seq)] * init_pop_size
        O.create_mtdna(init_pop_size)
        for m in O.mtdna_lst:
            #m.set_ttl(np.random.randint(1, 10))
            m.set_mtdna_size(16569)
            m.set_ttl(10)


        #O = msim.Mito(meta_data,mito_params)
        #attr_lst = [(wild_seq, mtdna_params)] * init_pop_size
        #O.create_mtdna(attr_lst)
        #msg = "Create oragnelle: %s" % (O.mito_id)
        #msim.write_debug(msg, 1)

        #for m in O.mtdna_lst:
        #    m.set_ttl(np.random.randint(1, 10))

        #O.display_params()

        results_lst.append(O.state)


        # Run simulation
        t1,t2 = O.run_sim(run_length)
        msg = "Run time: %s" % (t2 - t1)
        msim.write_debug(msg, 1)


        O.write_ts(ts_file, sim_id)
        O.write_gen(gen_file)


    for r in results_lst:
        print(r)

    with open(results_file, 'w') as fd:
        fd.writelines(json.dumps(results_lst))

    msg = "Results written to: %s" % (results_file)
    msim.write_debug(msg, 1)

