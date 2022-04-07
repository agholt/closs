#!/usr/bin/env python3

__author__ = "Alan Holt"
__copyright__ = "Copyright 2019"
__credits__ = ["Adrian Davies"]
__license__ = "Proprietary"
__version__ = "2.0.0"
__maintainer__ = "Alan Holt"
__email__ = "agholt@gmail.com"
__status__ = "Reduced GC"


import itertools
import datetime
import json
from io import StringIO
#import hashlib
import timeit
import csv
import gc

import random
import numpy as np
import cProfile

from mtDNA import *

DEBUG = int(3)


def proc_results(O):

    m = {'test' : 'hello'}
    wild_type = O.wild_type
    
    
    mut_lst = []
    t_start = int(O.species_dct[wild_type][0])
    t_end = int(len(O.species_dct[wild_type][2]))
    
    surv_time  = (t_start + t_end) - 1
    m['surv_time'] = int(surv_time)
    
    
    m["wild_start"] = int(O.species_dct[wild_type][1][t_start])
    m["wild"] = int(O.species_dct[wild_type][2][t_end - 1])
    
    mut_start_lst = []
    mut_end_lst = []
    k_lst = []
    
    for k in O.species_dct:
        
        if k != wild_type:
            if O.species_dct[k][0] == t_start:
                mut_start_lst.append(str(O.species_dct[k][1][0]))
    
            t1 = O.species_dct[k][0]
            t2 = len(O.species_dct[k][1])
            if (t1 + t2) == (t_end):
                mut_end_lst.append(O.species_dct[k][1][-1])
                k_lst.append(k)
                
    m['mut_start'] = " ".join(mut_start_lst)
    m['mutant_count'] = int(np.sum(mut_end_lst))

    mx = int(np.max(mut_end_lst))
    m['max_mutant'] = mx
    
    for k,s in zip(k_lst, mut_end_lst):
        if s == mx:
            biggest_mut = k
            break
            
    m['biggest_mut'] = str(k)
    for i in O.mtdna_lst:
        if i.type == biggest_mut:
            m['biggest_mut_size'] = int(i.sequence_len)
            break
    
    return m


def write_debug(msg, debug_level):
        
    if debug_level <= DEBUG:
        ts = str(datetime.datetime.now()).split('.')[0].replace(' ', '_')
        print("%s: %s" % (ts,msg))


class Common_methods():
    
    def debug_output(self, msg, debug_level):
        
        if self.debug >= debug_level:
            ts = str(datetime.datetime.now()).split('.')[0].replace(' ', '_')
            print("%s: %s" % (ts,msg))
            
        return None


class ATP():
        
    def init_atp(self, atp_level):
        self.atp = atp_level
        #self.atp_lst.append(self.atp)
        return self.atp

    def add_atp(self, atp_level):
        self.atp += atp_level
        return self.atp
    
    def del_atp(self, atp_level):

        consumed = min(self.atp, atp_level)
        self.atp -= atp_level
        self.atp = max(0,self.atp)
        #self.atp_lst.append(self.atp) # record atp level
        return consumed
    
    #def record_atp(self, atp_level, src, clone):
    #    
    #    self.atp_df['TIME'].append(self.t)
    #    self.atp_df['SRC'].append(src)
    #    self.atp_df['ATP'].append(atp_level)
    #    self.atp_df['CLONE'].append(clone)
        
        return None


class MITO_methods():
    
    def get_ID(self):
        self.ID_count += 1
        return self.ID_count
    
    def set_ID(self, count):
        self.ID_count = count + 1
        return None
    
    def set_species_db(self, mtseqs):
        self.mtseqs = mtseqs
        return None

        
    def add_mt (self, ID, param_dct):
        p = param_dct
        if not ID in self.mtids:
            p.update({'ID': ID})
            p.update({'tSTART': self.t})
            p.update({'tEND': 0})        
            self.mtids[ID] = p
    
    def set_mt_end_date(self, ID):
        if ID in self.mtids:
            self.mtids[ID]['tEND'] = self.t     

    def wipe_out(self):
        if self.wipe_out_count > 0:
            return
        if random.uniform(0, 1) < self.wipe_out_prob:
            msg = "wipe out: %s" % (self.t)
            self.debug_output(msg,0)
              
            for m in self.mtdna_lst:
                if random.uniform(0, 1) < self.wipe_out_rate:
                    m.ttl = 0
        self.wipe_out_count = self.wipe_out_hold 
            
    
    # Delete mtdna if TTL=0
    def delete_dead(self):
        i = 0
        for m in self.mtdna_lst:
            if m.ttl < 1:
                self.mtdna_lst.pop(i)
                self.set_mt_end_date(m.ID)

                m_key = self.get_key(m.species)
                self.mtdna_data[m_key][2][-1] -= 1
                

            i+=1

    def species_lst_str(self, s_lst):
        return "-".join([str(i) for i in s_lst])

    def species_str_lst(self, s_str):
        return s_str.split("-")


    def count_species(self):
        #self.count_species1()
        self.count_species3()
        return None

    def count_species1(self):

        species = [self.species_lst_str(m.species) for m in self.mtdna_lst]
        #species = ["-".join([str(i) for i in m.species]) for m in self.mtdna_lst]

        no_mtdna, no_wt, no_mt = 0,0,0
        pop_count = {}

        for s in species:
            no_mtdna += 1
            if s == "0":
                no_wt +=1
            else:
                no_mt +=1

            #s_ = "-".join([str(i) for i in s])
            if s in pop_count:
                pop_count[s] += 1
            else:
                pop_count[s] = 1

        self.count = no_mtdna
        self.count_lst.append(self.count)
        self.no_wild_type = no_wt
        self.no_mutants = no_mt

        self.current_no_species = len(set(species))
        self.no_species.append(self.current_no_species)

        for mtype in pop_count:

            count = pop_count[mtype]

            if mtype in self.species_dct:
                self.species_dct[mtype][2].append(count)
            else:    
                for m in self.mtdna_lst:
                     ms = "-".join(str(i) for i in m.species)
                     if ms == mtype:
                          coop_flag = m.cooperate_flag
                          break
                self.species_dct[mtype] = (self.t, coop_flag, [count])


        self.h[self.t] = (self.count,self.no_wild_type,self.no_mutants,self.current_no_species)



        return None


    # Count the numner of mtdna and the number of different species
    def count_species2(self):

        self.count = len(self.mtdna_lst)
        self.count_lst.append(self.count)

        
        #species = ["-".join([str(i) for i in m.species]) for m in self.mtdna_lst]
        species = [self.species_lst_str(m.species) for m in self.mtdna_lst]
        s_count = [(s, species.count(s)) for s in set(species)]
        pop_count = dict(s_count)

        self.no_wild_type = sum([s=="0" for s in species])
        self.no_mutants = self.count - self.no_wild_type


        self.current_no_species = len(s_count)
        self.no_species.append(self.current_no_species)

        coop_dct = {}
        for m in self.mtdna_lst:
            s = self.species_lst_str(m.species)
            if s not in coop_dct:
                coop_dct[s] = m.cooperate_flag

        
        for mtype in pop_count:
            count = pop_count[mtype]
            if mtype in self.species_dct:
                self.species_dct[mtype][2].append(count)
            else:    
                coop_flag = coop_dct[mtype]
                self.species_dct[mtype] = (self.t, coop_flag, [count])
                #print("debug 3:", coop_dct)
                #print("debug 4:", self.t, mtype, coop_flag, [count])

        self.h[self.t] = (self.count,self.no_wild_type,self.no_mutants,self.current_no_species)
                
        return None



    def count_species3(self):

        self.count = len(self.mtdna_lst)

        self.no_wild_type = 0
        self.no_mutants = 0
        self.current_no_species = 0 

        for k in self.mtdna_data:
            if k == "0":
                self.no_wild_type += self.mtdna_data[k][2][-1]
            else:
                self.no_mutants += self.mtdna_data[k][2][-1]

            t0 = self.mtdna_data[k][0]
            t1 = len(self.mtdna_data[k][2])
            if self.t == (t0+t1):
                self.current_no_species += 1


        
        self.no_species.append(self.current_no_species)

        self.h[self.t] = (self.count,self.no_wild_type,self.no_mutants,self.current_no_species)
            




    def perform_ufission(self):
        if self.ufission_frag_method == 'constant':
            self.ufission_const()
        if self.ufission_frag_method == 'poisson':
            self.ufission_poisson()
        if self.ufission_frag_method == 'percent':
            self.ufission_percent()

        return self.ufission_frag_method

    def ufission_const(self):

        
        if self.ufission_period == 0:
            return None

        if self.t % self.ufission_period == 0:
        #if random.uniform(0, 1) < self.ufission_bern_thresh:
            np.random.shuffle(self.mtdna_lst)
            frag = self.mtdna_lst[:self.ufission_const_frag_size]
            frag_size = len(frag)
            #msg = "ufission_const: %s %s" % (self.t, frag_size)
            #self.debug_output(msg,j)
            #write_debug(msg,4)

            no_wild = 0
            no_mut = 0
            for m in frag:
                if m.type == self.wild_type:
                    no_wild += 1
                else:
                    no_mut += 1
                # Here

 

            frag_health = float(no_mut)/float(frag_size) # high is bad

            if frag_health > self.ufission_health_thresh:
                new_lst = self.mtdna_lst[self.ufission_const_frag_size:]
                self.mtdna_lst = new_lst
                #msg = "ufission_const: %s (%s %s %s) fragment died" % (self.t, no_wild, no_mut, frag_size)
                #self.debug_output(msg,3)
                write_debug(msg,4)

                # decrement species count accordingly
                for m in frag:
                    m_key = self.get_key(m.species)
                    self.mtdna_data[m_key][2][-1] -= 1

           
        return None

    def ufission_poisson(self):

        self.poisson_flag = True

        #if self.t % self.ufission_period == 0:
        if random.uniform(0, 1) < self.ufission_bern_thresh:
            
            # Fission failure
            self.ufission_prob = np.minimum(0.95, self.ufission_fail_prob + (self.ufission_prob_per_mutant * self.no_mutants))
            if random.uniform(0, 1) < self.ufission_prob:
                #self.ff_rec.append((self.t,"ff_fail", "%d %f" % (self.no_mutants,self.ufission_prob)))
                return None  # ufission fails

            # Here
            np.random.shuffle(self.mtdna_lst)
            frag_size = np.random.poisson(self.ufission_poisson_mean_frag_size, size=1)[0]
            #a,b = self.ufission_frag_prop_limits
            #frag_size = int(random.uniform(float(a), float(b)) * self.count)
            #msg = "%s %s %s" % (a,b,frag_size)
            #self.debug_output(msg,1)

            frag = self.mtdna_lst[:frag_size]
            frag_size = len(frag)
        
            # fragment empty
            if frag_size == 0:
                return None


            # compute health of fragment 
            m_count = sum([m.species != [0] for m in frag])
            frag_health = float(m_count)/float(frag_size) # high is bad

            if frag_health > self.ufission_health_thresh:
                new_lst = self.mtdna_lst[frag_size:]
                self.mtdna_lst = new_lst
                #self.ff_rec.append((self.t, "fission", "%d %d %f" % (self.no_wild_type,frag_size,frag_health)))

                #msg = "%s (%s:%s) fragment died %d %d" % (self.t, self.no_wild_type, self.count, frag_size, m_count)
                #msg = "no fusion: [%5.4f] -- t: %s: w: %s m: %s -- frag size: %d no. muts: %d" % (self.t/float(96.0), self.count, self.no_wild_type, self.no_mutants, frag_size, m_count)

                #self.debug_output(msg,1)

        return None


    def ufission_percent(self):

        self.poisson_flag = True

        if random.uniform(0, 1) < self.ufission_bern_thresh:
            
            # Fission failure
            self.ufission_prob = np.minimum(0.95, self.ufission_fail_prob + (self.ufission_prob_per_mutant * self.no_mutants))
            if random.uniform(0, 1) < self.ufission_prob:
                return None  # ufission fails

            np.random.shuffle(self.mtdna_lst)
            a,b = self.ufission_frag_prop_limits
            frag_size = int(random.uniform(float(a), float(b)) * self.count)

            frag = self.mtdna_lst[:frag_size]
            frag_size = len(frag)
        
            # fragment empty
            if frag_size == 0:
                return None


            # compute health of fragment 
            m_count = sum([m.species != [0] for m in frag])
            frag_health = float(m_count)/float(frag_size) # high is bad

            if frag_health > self.ufission_health_thresh:
                new_lst = self.mtdna_lst[frag_size:]
                self.mtdna_lst = new_lst

                # decrement species count accordingly
                for m in frag:
                    m_key = self.get_key(m.species)
                    self.mtdna_data[m_key][2][-1] -= 1

        return None



        
    def save_living_mtdnas(self):
        
        self.time_to_fission_thresh = self.t
        for m in self.mtdna_lst:
            m_rec = {}
            m_rec['ID'] = m.ID
            m_rec['PID'] = m.PID
            m_rec['TYPE'] = "-".join(str(m) for m in m.species)
            m_rec['TTL'] = m.ttl
               
            m_rec['ISCLONE'] = m.is_clone 
        
            m_rec['BUSY'] = m.busy 
            m_rec['GEN'] = m.generation
            m_rec['ATP_FLAG'] = m.ATP_FLAG 
        
            self.living_mtdnas.append(m_rec)
            
        return None
    
    def get_living_mtdnas(self):
        
        self.time_to_fission_thresh = self.t
        living_mtdnas = []
        for m in self.mtdna_lst:
            m_rec = {}
            m_rec['ID'] = m.ID
            m_rec['PID'] = m.PID
            #m_rec['TYPE'] = m.type
            m_rec['TYPE'] = "-".join(str(m) for m in m.species)
            m_rec['TTL'] = m.ttl
               
            m_rec['ISCLONE'] = m.is_clone 
        
            m_rec['BUSY'] = m.busy 
            m_rec['GEN'] = m.generation
            m_rec['ATP_FLAG'] = m.ATP_FLAG 
        
            living_mtdnas.append(m_rec)
            
        return living_mtdnas

    def _get_exp_params(self):

        mutate_model_params = (0)
        if self.mutate_model == "LIN":
            mutate_model_params = self.lin_params
        if self.mutate_model == "EXP":
            mutate_model_params = self.exp_params
        if self.mutate_model == "EXP2":
            mutate_model_params = self.exp2_params
        if self.mutate_model == "PW":
            mutate_model_params = self.pw_params

        params = {
           'capacity': self.max_population,
           'ufission_health_thresh': self.ufission_health_thresh,
           'ufission_poisson_mean_frag_size': self.ufission_poisson_mean_frag_size,
           'ufission_const_daily_rate': self.ufission_const_daily_rate,
           'ufission_fail_prob': self.ufission_fail_prob,
           'ufission_prob_per_mutant': self.ufission_prob_per_mutant,
           'mutate_prob_per_mutant': self.mutate_prob_per_mutant,
           'mutate_model': self.mutate_model,
           'mutate_model_params': mutate_model_params 
        }

        return params

        
    def simulation_output(self, t,u0):

        # Periodic printout            
        if t % (24 * self.sample_per_hour * self.day_output_rate) == 0:

            u1 = timeit.default_timer()
            if t == 0:

                print("_START SIMULATION_")
                print("run_length: %s" % (self.run_length))
                print(self.experimental_params)

                return u1

            delay = u1 - u0
            day = t/float(24 * self.sample_per_hour)
            #age_yr = float(self.age)/float(day * 365)
            age_yr = float(self.age_offset + t)/float(24 * self.sample_per_hour * 365)

            msg = "days: %d [%.2f] t: %d  w: %d m: %d s :%d -- age %.2f | mprob %.6f | fis fail prob %.2f" %  \
                 (int(day),  
                  delay, 
                  self.count, 
                  self.no_wild_type, 
                  self.no_mutants, 
                  self.current_no_species, 
                  age_yr, 
                  self.mutate_prob, 
                  self.ufission_prob)

            self.debug_output(msg,1)

            #msg = "alternate count - t: %d w: %d m: %d" %  (self.total_count, self.wt_count, self.mt_count)
            #self.debug_output(msg,1)

            return u1

        return u0

    # Write generation data
    def write_gen(self, gen_file):
        #gen_file = "mito-gen_%s.json" % (ts)
        #gen_file = "mito-gen.json"
    
        mtids_lst = [self.mtids[k] for k in self.mtids]
        df = {"META": {}, "DATA": {}}
        df['META'] = {'wild_type': self.wild_type, 'run_length': self.run_length, 'samples_per_day': self.sample_per_hour * 24}
        df['DATA'] = self.mtids
  
        with open(gen_file, 'w') as fd:
            fd.writelines(json.dumps(df))
        msg = "Generation results written to: %s" % (gen_file)
        write_debug(msg, 2)
    
    
    # Write timeseries data
    def write_ts(self, ts_file, sim_id):
        
        df = {"META": {}, "DATA": {}}
        df['META'] = {'wild_type': self.wild_type, 
                      'run_length': self.run_length, 
                      'age_start': self.age_start, 
                      #'time_to_fission_thresh': self.time_to_fission_thresh, 
                      'sample_per_day': 24 * self.sample_per_hour}
        df['DATA'] = self.species_dct
        df['H'] = self.h
        df['TOTAL'] = self.count_lst
        df['NO_SPECIES'] = self.no_species
        #df['ATP'] = self.atp_lst
        df['SIM_ID'] = sim_id

        df['PARAMS'] = self.experimental_params

        with open(ts_file, 'w') as fd:
            fd.writelines(json.dumps(df))
        msg = "Timeseries results written to: %s" % (ts_file)
        write_debug(msg,1)

    def ts_save(self, ts_file):
        df = {}
        #print(self.mtdna_data)
        df['DATA'] = self.mtdna_data
        df['META'] = {'wild_type': self.wild_type, 
                      'run_length': self.run_length, 
                      'age_start': self.age_start, 
                      'sample_per_day': 24 * self.sample_per_hour}
        df['PARAMS'] = self.experimental_params
        df['H'] = self.h

        with open(ts_file, 'w') as fd:
            fd.writelines(json.dumps(df))
        msg = "Timeseries results written to: %s" % (ts_file)
        write_debug(msg,1)

    def write_state_file(self, state_file, meta_params, mito_params, mtdna_params, sim_id):

        max_mid = np.max([m['ID'] for m in self.living_mtdnas])

        meta_params['run_length'] = self.run_length, 
        #meta_params['time_to_fission_thresh'] = self.time_to_fission_thresh, 

        mito_state = {}
        #mito_state['MAX_MID'] = max_mid
        mito_state['META_PARAMS'] = meta_params
        mito_state['MITO_PARAMS'] = mito_params
        mito_state['MTDNA_PARAMS'] = mtdna_params
        mito_state['MTDNA_LST'] = self.living_mtdnas
        mito_state['SPECIES_DB'] = self.mtseqs
        mito_state['WILD_TYPE'] = self.wild_type
        mito_state['SIM_ID'] = sim_id

        with open(state_file, 'w') as fd:
            fd.writelines(json.dumps(mito_state))

        msg = "Mito state written to: %s" % (state_file)
        write_debug(msg,1)
        return None


    # Set parameters methods
    def _get_param2(self, params, param_key, default_val):
        
        if param_key in params:
            new_param = params[param_key]
        else:
            new_param = default_val
                                    
        return new_param

    def _set_admin_params(self):

        p = self.params['ADMIN_PARAMS']


        self.run_id = self._get_param2(p, 'RUN_ID', 1)
        self.day_output_rate = self._get_param2(p, 'DAY_OUTPUT_RATE', 60)
        self.VERBOSE = self._get_param2(p, 'VERBOSE', True)
        self.debug = self._get_param2(p,'DEBUG_LEVEL', 0)
        self.gen_logging_on = self._get_param2(p,'GEN_LOGGING_ON', False)
        self.log_to_file = self._get_param2(p,'LOG_TO_FILE', False)
        self.log_file = self._get_param2(p,'LOG_FILE', "msim.log")
        self.log_only_to_file = self._get_param2(p,'LOG_ONLY_TO_FILE', False)
        self.terminate_on_extinction = self._get_param2(p,'TERMINATE_ON_EXTINCTION', True)

        return None

    def _set_init_params(self):

        p = self.params['INIT_PARAMS']

        self.max_population = self._get_param2(p, 'max_population', None)
        self.age_start = self._get_param2(p, 'age_start', None)
        self.wild_type_thresh = self._get_param2(p, 'WILD_TYPE_THRESH', None)

        return None


    def _set_mtdna_params(self):

        p = self.params['MTDNA_PARAMS']


        # Busy and TTL paramaters
        self.max_ttl = self._get_param2(p, 'MAX_TTL', None)
        self.max_repl_time = self._get_param2(p,'MAX_REPL_TIME', False)
        self.constant_repl_time = self._get_param2(p, 'CONSTANT_REPL_TIME', False)

        # Event probabilities
        self.damage_prob = self._get_param2(p, 'DAMAGE_PROB', None)

        clone_probs = self._get_param2(p, 'CLONE_PROBS', (0.0,0.0,0.0))
        self.p1, self.p2, self.p3 = clone_probs

        # Clone probability
        # replS, either used as static clone probability or as part of Verhulst pop. eqn.
        self.replS = self._get_param2(p, 'CLONE_replS', None)
        self.clone_scale_flag = self._get_param2(p, 'CLONE_SCALE_FLAG', False)

        self.clone_prob_method = self._get_param2(p, 'CLONE_PROB_METHOD', None)

        self.no_clone_without_atp = self._get_param2(p, 'NO_CLONE_WITHOUT_ATP', False)

        self.defect_on = self._get_param2(p, 'DEFECT_ON', 0)

        self.wipe_out_prob = self._get_param2(p, 'WIPE_OUT_PROB', 0)
        self.wipe_out_rate = self._get_param2(p, 'WIPE_OUT_RATE', 0)
        self.wipe_out_hold = self._get_param2(p, 'WIPE_OUT_HOLD', 0)

        # Binge drinking
        self.binge_cycle = self._get_param2(p, 'BINGE_CYCLE', 0)
        self.binge_attr = self._get_param2(p, 'BINGE_ATTR', 0.5)



        return None
        

    # Set fission parameters
    def _set_fission_params(self):

        p = self.params['FISSION_PARAMS']

        self.ufission_frag_method = self._get_param2(p, 'ufission_frag_method', None)
        self.ufission_rate_method = self._get_param2(p, 'ufission_rate_method', None)

        self.ufission_health_thresh = self._get_param2(p, 'ufission_health_thresh', None)
        self.ufission_const_frag_size = self._get_param2(p, 'ufission_const_frag_size', None)
        self.ufission_poisson_mean_frag_size = self._get_param2(p, 'ufission_poisson_mean_frag_size', None)
        self.ufission_frag_prop_limits = self._get_param2(p, 'ufission_frag_prop_limits', None)
        self.ufission_const_daily_rate = self._get_param2(p, 'ufission_const_daily_rate', None)

        self.ufission_fail_prob = self._get_param2(p, 'ufission_fail_prob', None)
        self.ufission_prob_per_mutant = self._get_param2(p, 'ufission_prob_per_mutant', None)

        if self.ufission_const_daily_rate == 0:
            self.ufission_period = 0 
        else:
            self.ufission_period = int((24.0 * self.sample_per_hour)/float(self.ufission_const_daily_rate) + 0.5)

        return None


    def _set_freerad_params(self):

        p = self.params['FREERAD_PARAMS']

        self.mutate_prob_base = self._get_param2(p, 'mutate_prob_base', 1e-6)
        self.non_replicable_prop = self._get_param2(p, 'non_replicable_prop', 0.0)
        self.defection_prob = self._get_param2(p, 'defection_prob', 0.0)
        self.mutate_prop_flag = self._get_param2(p, 'mutate_prop_flag', False)
        self.mutate_prob_per_mutant = self._get_param2(p, 'mutate_prob_per_mutant', 1e-6) # /float(self.sample_per_hour)

        self.const_model = self._get_param2(p, 'const_model', 1e-6)
        self.mutate_model = self._get_param2(p, 'mutate_model', None)
        self.lin_params = self._get_param2(p, 'lin_model', (None,None))
        self.exp_params = self._get_param2(p, 'exp_model', (None,None,None))
        self.pw_params = self._get_param2(p, 'pw_model', (None,None,None,None,None))

        return None


    def _set_atp_params(self):

        p = self.params['ATP_PARAMS']

        self.consume_rate = self._get_param2(p,'atp_consume_rate', None)
        self.atp_seq = self._get_param2(p,'ATP_SEQ', None)
        self.produce_rate = self._get_param2(p,'atp_produce_rate', None)       # mtDNA ATP production
        self.repl_consume_rate = self._get_param2(p,'repl_atp_consume_rate', None)
        self.per_mtdna_atp_consumuption = self._get_param2(p,'per_mtdna_atp_consumuption', 0)
        self.atp_production_always_on  = self._get_param2(p,'ATP_PRODUCTION_ALWAYS_ON', True)
        self.allow_mutant_atp_producers  = self._get_param2(p,'ALLOW_MUTANT_ATP_PRODUCERS', False)
 
        return None


    def _set_sim_time_params(self):

        p = self.meta_data['sim_time_params']
        self.sample_per_hour = self._get_param2(p, 'sample_per_hour', None)

        return None




class Mito(Common_methods, MITO_methods, ATP):
        
    def __init__(self, meta_data, params):

        ts = str(datetime.datetime.now()).split('.')[0].replace(' ', '_')
        self.mito_id = str(hashlib.sha224(ts.encode('utf-8')).hexdigest())

        self.params = params
        self.meta_data = meta_data

        #self._set_admin_params2() # Legacy params
        self._set_admin_params()

        self._set_mtdna_params()

        self._set_sim_time_params()
        self._set_init_params()
        self._set_fission_params()
        self._set_freerad_params()
        self._set_atp_params()

        self._init_data()

        self.non_rep_total = 0
        self.mtdna_total = 0
        self.clone_events = 0

        self.mutate_count = 0



        self.h = {}

        if self.VERBOSE:
            print("Author: %s"  % (__author__))
            print("Credits: %s" % (__credits__))
            print("Version: %s" % (__version__))
            print("Status: %s"   % (__status__))
            print("Copyright: %s" % (__copyright__))

            print(self.params['FREERAD_PARAMS'])
        
    
    def _init_data(self):
         
        self.ID_count = 0
        self.CLONE_FLAG = 1

        self.wipe_out_count = self.wipe_out_hold

        self.gc_rate = 1
            
        self.t = 0
        self.age_offset = int(self.age_start * 24 * self.sample_per_hour * 365)

        self.mutate_prob = self.mutate_prob_base # (self.m * self.age/(24 * self.sample_per_hour * 365.0)) + self.c 
        
        self.mtda_lst = []
        self.state = {}
        
        self.count = int(0)
        self.no_wild_type = int(0) 
        self.no_mutants = int(0) 
        self.mutant_count = 0
        self.count_lst = []

        self.total_count = int(0)
        self.wt_count = int(0)
        self.mt_count = int(0)
        
        self.mtseqs = {}
        self.mtids = {}
        
        self.species_lst = []
        self.species_dct = {}

        #self.spec_lst = [[0]]
        
        self.no_species = []
        
        self.atp = int(0)
        self.atp_lst = []        
        self.atp_producer_lst = []
        
        self.m_state = []
        
        self.living_mtdnas = []
        self.fission_thresh_flag = False
 
        self.ufission_prob = float(0.0)

        self.ff_rec = []


        self.mutant_prob = 0.0


        self.ufission_bern_thresh = self.ufission_const_daily_rate/(24.0 * self.sample_per_hour)

        self.defector_lst = []


        self.wild_type = 0 
        self.wild_len = 16569 


        self.defector_start_time = 0

        return None

    def inc_counts(self, m):
        self.total_count  += 1
        if m.type == self.wild_type:
            self.wt_count += 1
        else:
            self.mt_count += 1

    def dec_counts(self, m):
        self.total_count  -= 1
        if m.type == self.wild_type:
            self.wt_count -= 1
        else:
            self.mt_count -= 1

    def get_key(self, species_lst):
       return "-".join([str(i) for i in species_lst])


    
    def create_mtdna(self, n):

        self.mtdna_lst = []
        for i in range(n):
            m = mtDNA(self)
            m.ID = self.get_ID() 
            m.PID = 0
            m.busy = 0
            m.generation = 0
            m.species = [0]
            m.set_spec_level()
            #m.mtdna_size = 16569
            m.set_mtdna_size(16569)
            m.update_env()
            self.mtdna_lst.append(m)

            #self.inc_counts(m)
        
        wild_spec = [0]    
        wild_spec_key = self.get_key(wild_spec)
        self.spec_lst = [wild_spec]
        self.spec_str_lst = [wild_spec_key]

        # mtdna_data created
        self.mtdna_data = { wild_spec_key: (self.t, True, [n]) }


        self.count_species()

        return None

    def add_mtdna(self,seq):

        m = mtDNA(seq, self)
        m.ID = self.get_ID() 
        m.PID = 0
        m.busy = 0
        m.generation = 0
        m.update_env()
        self.mtdna_lst.append(m)

        self.count_species()
        return None
        
    
    #def save_state2(self, m_dct):
    #    self.m_state.append(m_dct)
        
    def save_state(self, reason):

        self.state['DESCRIPTION'] = reason
        self.state['TIME'] = self.t
        self.state['RUN_LENGTH'] = self.run_length
        self.state['DEFECT_START_TIME1'] = self.defector_start_time
        #self.state['DEFECT_START_TIME'] = self.mtdna_data['0-0'][0]
        #self.state['ATP'] = self.atp_lst[-1]
        #self.state['NO_MTDNA'] = self.count_lst[-1]
        self.state['NO_MTDNA'] = self.count
        self.state['NO_WILD_TYPE'] = self.no_wild_type
        self.state['NO_MUTANTS'] = self.no_mutants
        self.state['NO_MUTANT_SPECIES'] = self.current_no_species


        self.state['TOTAL_SPECIES'] = len(set(self.mtdna_data.keys()))
        self.state['ufission_frag_method'] = self.ufission_frag_method
        self.state['ufission_period'] = self.ufission_period
        self.state['ufission_health_thresh'] = self.ufission_health_thresh


        self.state['FF_fail_per_mutant'] = self.ufission_prob_per_mutant
        self.state['MUTATION_PROB'] = self.const_model
        self.state['MUTATION_PROB_per_mutant'] = self.mutate_prob_per_mutant



        #if self.ufission_frag_method == 'constant':
        #    self.state['ufission_frag_size'] = self.ufission_const_frag_size
        #elif self.ufission_frag_method == 'poisson':
        #    self.state['ufission_frag_size'] = self.ufission_poisson_mean_frag_size
        #elif self.ufission_frag_method == 'percent':
        #    self.state['ufission_frag_size'] = self.ufission_frag_prop_limits
        #else:
        #    self.state['ufission_frag_size'] = None

        #self.state['EXPR_PARAMS'] = self.experimental_params

        self.state['DATETIME'] = str(datetime.datetime.now()).split('.')[0].replace(' ', '_')

        return None

    def binge(self,t, bing_cyc):
        if t == 0:
            return 0
        if (t % bing_cyc) == 0:
            L = len(self.mtdna_lst)
            np.random.shuffle(self.mtdna_lst)
            L2 = int(float(L) * self.binge_attr)
            self.mtdna_lst[L2:] = []
            
            print("[%s] binge cycle. No. mtdna %d / %d" % (t, L, L2))

    
    # Run the simaultion
    def run_sim(self, run_length):

        self.run_length = run_length

        self.experimental_params = self._get_exp_params()

        t1 = u0 = timeit.default_timer()


        bing_cyc = self.binge_cycle * 24 * self.sample_per_hour
        print("BINGE CYCLE:",  self.binge_cycle, bing_cyc)
        
        for t in range(self.run_length):
            
            self.t += 1

            child_lst = []

            # Print simulation output
            u0 = self.simulation_output(t,u0)

        

            # Initialise mtdna counts for next iteration
            for k in self.mtdna_data:
                count = self.mtdna_data[k][2][-1]
                if count > 0:
                    self.mtdna_data[k][2].append(count)

            self.wipe_out_count +- 1
            self.wipe_out()

            # binge drinking
            if self.binge_cycle > 0 and t > 0:
                if (t % bing_cyc) == 0:
                    print("[%s] binge drink" % (t))      
                    for m in self.mtdna_lst:
                        if random.uniform(0, 1) < self.binge_attr:
                            m.ttl = 0

                

            # Purge dead mtdna
            self.delete_dead()

            consumed = self.del_atp(self.consume_rate)  

            # This is a bit of a hack:
            self.mtdna_lst[0].calc_mutation_probability()

            # Perform microfission/fusion
            self.perform_ufission()
            
            #wild_type_count = mutant_count = 0

            for m in self.mtdna_lst:

                # Turn ATP production on/off
                #if consumed < self.consume_rate:
                #    m.set_atp(True)
                #else:
                #    m.set_atp(False)

                # Age and reduce process busy period (if in busy period)    

                m.calc_status()

        
                # Clone
                if consumed < self.consume_rate or (not m.cooperate_flag):
                    child = m.clone()
                    if child != None:
                        self.mtdna_lst.append(child)
                        k = self.get_key(child.species)

                        if k in self.mtdna_data:
                            self.mtdna_data[k][2][-1] += 1
                        else:
                            self.mtdna_data[k] = (self.t, child.cooperate_flag, [1])

                        #self.inc_counts(child)
 
                    
            # count species and other house keeping
            self.count_species()
 
            if t % (self.gc_rate * 364 * 24 * self.sample_per_hour) == 0:
                gc.collect()


            # Cell dies 
            if self.no_wild_type < self.wild_type_thresh:
                self.save_state("wild type threshold reached: %d" % (self.no_wild_type))
                self.time_to_wild_thresh = self.t

                if self.terminate_on_extinction:
                    t2 = timeit.default_timer()
                    msg = "wild type threshold reached at time %s no. wild type: %d" % (self.t, self.no_wild_type)
                    write_debug(msg,1)
                    return t1,t2
            
        # End of the simulation has been reached
        self.save_state("end of simulation")
        self.save_living_mtdnas()
        self.time_to_fission_thresh = self.t
        t2 = timeit.default_timer()
        msg = "end of simulation %s" % (self.t)
        write_debug(msg,1)

        return (t1,t2)





