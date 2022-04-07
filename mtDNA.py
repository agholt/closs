#!/usr/bin/env python3

__author__ = "Alan Holt"
__copyright__ = "Copyright 2019"
__credits__ = ["Adrian Davies"]
__license__ = "Proprietary"
__version__ = "1.6.0"
__maintainer__ = "Alan Holt"
__email__ = "agholt@gmail.com"
__status__ = "Defection"


import itertools
import datetime
import json
from io import StringIO
import hashlib
import timeit
import csv

import random
import numpy as np

DEBUG = 3

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



class mtDNA_methods():
    
    
    # Return the mtdna's ID
    def get_id(self):
        return self.ID
    
    # Set the mtdna's parent ID
    def set_parent(self, pid):
        self.PID = pid
        return self.PID
    
    # Return the parent ID of the mtdna
    def get_parent(self):
        return self.PID
    
    # start replication timer (called from clone method)
    def set_busy(self, busy):
        self.busy = busy
        return None
    
    def set_as_clone(self):
        self.is_clone = True
        return None

    def set_spec_level(self):
        self.species_level = len(self.species)

    def set_mtdna_size(self, siz):
        self.mtdna_size = siz
        self.mtdna_scale = float(self.mtdna_size)/float(self.env.wild_len)
        self.replication_prob = self.env.replS * self.mtdna_scale

        self.busy_scale = self.mtdna_scale
        self.max_busy = max(1, int(self.mtdna_scale * self.env.max_repl_time))
 
   
    # damage has occured, decrement TTL  
    def dec_ttl(self):
        self.ttl -= 1
        self.ttl = max(0, self.ttl)
        return None
    
    # Decrement busy period
    def dec_busy(self):
        self.busy -= 1
        self.busy = max(0, self.busy)
        if self.busy == 0:
            self.is_clone = False
        return None

    # Set time to live
    def set_ttl(self, ttl):
        self.ttl = ttl
        return None
        
    # Computer next generation
    def inc_generation(self,last_gen):
        self.generation = last_gen + 1
        return None
    
    def rnd_event(self, lbound,ubound, thresh):
        return np.random.randint(lbound, ubound) < thresh
    
    # Computes the mtdna'a busy period
    #def compute_max_busy(self):
    #    
    #    self.busy_scale = float(self.mtdna_size)/float(self.env.wild_len)
    #    self.max_busy = max(1, int(self.busy_scale * self.env.max_repl_time))

    #    return None
    
    # Update environment
    def update_env(self):
        
        if self.env.gen_logging_on:
            mparams = {
                'PID': self.PID,
                'TYPE': self.type,
                 'SEQ_LEN': self.mtdna_size,
                 'ATP_PRODUCER': self.atp_producer,
                 'GEN': self.generation
            }
            self.env.add_mt(self.ID, mparams)
        
        return None



class mtDNA(Common_methods, mtDNA_methods):
    
    def __init__(self, env):

        #self.params = params
        self.env = env
        self.env.mtdna_total += 1
        
        
        self.ID = self.env.get_ID() 
        self.PID = 0
        self.busy = 0
        self.generation = 0

        self.species = [0]
        self.spec_level = 0
        #self.set_mtdna_size = 16569
        
        self.ttl = self.env.max_ttl
        self.replicable_flag = True
        self.cooperate_flag = True


        self.replicable_mutate_prob = 0.0
        self.non_replicable_mutate_prob = 0.0


        #self.mutant_count = 0
        
        self.is_clone = False
        self.ATP_FLAG = True
        self.no_clone_without_atp = False
        
        #self.set_mtseq(mtseq) # add mitochondrial sequence to database
        
        #self.compute_max_busy()
        
        self.can_produce_apt()
        #if self.atp_producer:
        #    self.env.atp_producer_lst.append((self.ID,self.type))

        #self.mtScale = float(self.sequence_len)/float(self.env.wild_len)
        #self.replication_prob = self.env.replS * self.mtScale
        
    
    def can_produce_apt(self):
        
        self.atp_producer = False # initalise atp_producer to false
                                             
        # All wild types are ATP producers
        if self.species == [0]:
            self.atp_producer = True
            return None
        
        return None
        
    def not_replicable(self):
        self.replicable_flag = False
        

    # Set ATP generation flag (True of False)
    def set_atp(self, state):
        self.ATP_FLAG = state
        return None
    
    # Produce/consume ATP
    def produce_atp(self,atp_level):
        # If atp_production_always_on flag set, then override self.ATP_FLAG
        if self.atp_production_always_on:
            atp_flag = True
        else:
            atp_flag = self.ATP_FLAG
        
        if self.atp_producer and atp_flag:
            self.env.add_atp(atp_level)

        return None

   
    
    def consume_atp(self, atp_level):
        consumed = self.env.del_atp(atp_level)
        return consumed

    # Calculate mutation probability
    def calc_mutation_probability(self):

        ta = float(self.env.t + self.env.age_offset)/float(24 * self.env.sample_per_hour * 365)


        if self.env.mutate_model == 'CONST': 
             base_mprob = self.env.const_model

        if self.env.mutate_model == 'EXP2': # Not sure about this model...?
            a,b,c = self.env.exp2_params
            base_mprob = (ta * ta) * a * np.exp(b * ta + c) 
        if self.env.mutate_model == 'EXP':
            a,b,c,k = self.env.exp_params
            base_mprob = a * np.exp(b * ta + c) + k
        if self.env.mutate_model == 'LIN':
            m,c = self.env.lin_params
            base_mprob = m * ta + c 
        if self.env.mutate_model == 'PW':
            m1,c1,m2,c2,thresh = self.env.pw_params
            if ta < thresh:
                base_mprob = m1 * ta + c1 
            else:
                base_mprob = m2 * (ta-thresh) + c2 + (m1 * ta + c1)


        mprob = np.minimum(0.99, base_mprob + (self.env.mutate_prob_per_mutant * self.env.no_mutants))
        self.env.mutate_prob = mprob

        if self.env.mutate_prop_flag == True:
            #scale = float(self.sequence_len)/float(self.env.wild_len)
            mprob = self.mtdna_scale * mprob

        self.replicable_mutate_prob = mprob
        self.non_replicable_mutate_prob = self.replicable_mutate_prob / (1.0 - self.env.non_replicable_prop) 

        self.env.mutate_prob = self.replicable_mutate_prob 


        return self.replicable_mutate_prob, self.non_replicable_mutate_prob 


    
    # Calculate the status of the mtdna
    def calc_status(self):
        
        self.dec_busy() # calc busy time
        
       
        if self.species == [0]:
            self.env.add_atp(self.env.produce_rate)

        consumed = self.env.del_atp(self.env.per_mtdna_atp_consumuption)

        #if consumed < self.per_mtdna_atp_consumuption:
        #    self.dec_ttl()
        
        if random.uniform(0, 1) < self.env.damage_prob:
            self.dec_ttl()

        # Calculate mutation probability
        #self.calc_mutation_probability()

        return None

    def mutate(self):

        slen = len(self.species)
        slist = [s for s in self.env.spec_lst if len(s) == slen+1]
        print("debug1:", slist)

        next_spec = 0
        if slist != []:
            next_spec = np.max([s[-1] for s in slist]) + 1
            #print("debug2:", next_spec)
            
        new_spec = self.species.copy()
        new_spec.append(next_spec)            
        #print("debug3:", new_spec)

        #new_spec_flag = False
        #for i in range(10000000):
        #    new_spec = self.species.copy()
        #    new_spec.append(i)            
        #    if new_spec not in slist:
        #        break


        self.env.spec_lst.append(new_spec)
        new_size = int(self.mtdna_size/2.0)
        msg = "mutant: (%s) %s %s %s -- p: %5.4f" % (self.env.t, self.species, new_spec, new_size, self.replicable_mutate_prob)
        write_debug(msg,2)

        #print("species list:", self.env.spec_lst)

        return new_spec, new_size


    # Mutate
    def mutate2(self, mtseq):
        
        # Take a fragment of the mitochondrial sequence
        frag_len = 3 * int(self.no_codons/2.0)

        if frag_len > 1000:
            offset = np.random.randint(0, self.sequence_len-1)
            rotated_seq = mtseq[offset:] + mtseq[:offset] 
            new_seq = rotated_seq[frag_len:]
            msg = "mutate t: %s frag_len: %s" % (self.env.t, frag_len)
            write_debug(msg,4)
        else:
            new_seq = mtseq
        
        return new_seq

    
    # Compute clone probability as a function of population size
    def get_clone_prob(self):
        if self.env.clone_scale_flag:
            return self.replication_prob
        else:
            return self.env.replS 
        
    # clone mtdna
    def clone(self):

        self.env.clone_events += 1
        if not self.replicable_flag: # If not a repicable mutant return without cloning
            self.env.non_rep_total += 1
            #print("not replicable %s %s" % (self.ID, self.type))
 
            return None

        if self.env.count >= self.env.max_population:
            return None
        
        if self.ttl > 0 and self.busy == 0:
            
            clone_prob = self.get_clone_prob()
            
            if random.uniform(0, 1) < self.get_clone_prob(): # and self.env.count < self.env.max_population:
                
                self.busy = self.max_busy
                
                replicable_flag = True  # it's a replicable mutant, set flag
                defect_flag = False

                #new_seq = self.env.mtseqs[self.type]

                # Ensure no mutation beyond 2048 size
                if self.mtdna_size < 2000:
                    mut_thresh = 1.0
                else:
                    mut_thresh = random.uniform(0, 1) 

                # Calculate mutation probability
                self.calc_mutation_probability()




                new_spec = None
                if mut_thresh < self.replicable_mutate_prob:
                    #new_seq = self.mutate(new_seq)           
                    new_spec, new_size = self.mutate()           
                    replicable_flag = True  

                m_clone = mtDNA(self.env)      # clone

                m_clone.set_busy(self.max_busy) # set clone as busy
                m_clone.set_as_clone()
                m_clone.set_parent(self.ID)                   # set parent
                m_clone.inc_generation(self.generation)       # Put into next generation

                # Set inheritted traits:
                #if m_clone.type != self.type:

                m_clone.cooperate_flag = self.cooperate_flag

                # This is a mutant
                if new_spec != None:
                    m_clone.species = new_spec
                    m_clone.set_mtdna_size(new_size)
                    #self.env.spec_lst.append(m_clone.species)
                    #msg = "mutant: %s %s %s" % (m_clone.ID, m_clone.PID, m_clone.species)
                    #write_debug(msg,1)

                    if random.uniform(0, 1) < self.env.defection_prob:
                        m_clone.cooperate_flag = False
                        #self.env.defector_lst.append((m_clone.species,m_clone.ID))
                        msg = "defection: %s %s %s" % (m_clone.ID, m_clone.PID, m_clone.species)
                        write_debug(msg,1)

                # not a mutant
                else:
                    m_clone.species = self.species
                    m_clone.set_mtdna_size(self.mtdna_size)

                m_clone.set_spec_level()

                if m_clone.species not in self.env.species_lst:
                    self.env.species_lst.append(m_clone.species)


                # m_clone.update_env()

                
                return m_clone
            
        return None   


