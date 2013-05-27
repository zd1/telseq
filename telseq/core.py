'''
-----------------------------------------------
Copyright 2013 Wellcome Trust Sanger Institute
Written by Zhihao Ding (zd1@sanger.ac.uk)
Released under the GPL
-----------------------------------------------

Created on 26 May 2013

'''

import os
import sys
import re
import subprocess
import glob
import util
import pdb
    
__all__ = ['TelomereLength']

class SeqPattern():
        
    telomere_pattern = "TTAGGG"
    nullvalsymbol="NA"
    gc_interval=2
    
    def __init__(self, outdir, name, bamfile=None):
        
        self.patterns = [self.telomere_pattern]
        self.patterns.extend(self._tel_perm_letter())
        self.patterns = set(self.patterns)
        
        self.bamfile=bamfile
        self.bamproperty={
            "name":name,
            "bam":bamfile,
            "tol":0, # total reads 
            "map":0, # mapped reads
            "dup":0, # duplicates
            "gc":{}, # reads at in gc range
            "ptn":{}, # telomeric read counts for reads with different amount of patterns
        }
    
        for ptn in self.patterns:
            self.bamproperty['ptn'][ptn]={}
        
        self.outdir=outdir
        self.name=name
        
    def bamparser(self):
        
        # only extract flag, chromsome name, position, and the actual read
        samfields = "2,3,4,10"        
#         p1 = subprocess.Popen(["samtools","view", self.bamfile], stdout=subprocess.PIPE )
        p1 = subprocess.Popen(["cat", self.bamfile], stdout=subprocess.PIPE )
        p2= subprocess.Popen(["cut","-f%s"%samfields],stdin=p1.stdout, stdout=subprocess.PIPE )
        
        while (p2.poll() == None):
            if p2.stdout == None: break
            for line in p2.stdout:
                line = line.strip()
                if not line: break
                (flag, chrm, pos, read) = line.split("\t")
                self.bamproperty['tol'] += 1
                
                setflags = util.flagresolve(flag)
                if not setflags.has_key(0x4): # 0x4 for unmapped
                    self.bamproperty['map'] +=1
                if setflags.has_key(0x400): # 0x400 for pcr duplicates
                    self.bamproperty['dup'] +=1
                    
                gcfraction = util.gc_fraction(read)
                gcbin=int(gcfraction*100/self.gc_interval) # count read by its GC to every 2% GC bin
                
                if not self.bamproperty['gc'].has_key(gcbin):
                    self.bamproperty['gc'][gcbin]=1
                else:
                    self.bamproperty['gc'][gcbin]+=1
                
                for ptn in self.patterns:
                    target_count = self._countpattern(ptn, read)
                    if not self.bamproperty['ptn'][ptn].has_key(target_count):
                        self.bamproperty['ptn'][ptn][target_count] =1
                    else:
                        self.bamproperty['ptn'][ptn][target_count] +=1
                    
        util.cdm(self.bamproperty,"%s/%s.bamscan.pickle"%(self.outdir,self.name))
        
    @classmethod    
    def integrate(self, outdir):
        
        data={}
        dummy_name=0
        resultfiles = glob.glob("%s/*.bamscan.pickle"%outdir)
        if len(resultfiles) == 0:
            print >> sys.stderr, "No result files found in \n %s"%outdir
            sys.exit(1)
            
        for pkl in resultfiles:
            result = util.cl(pkl)
#             data[result['name']]=result
            data[dummy_name]=result
            dummy_name +=1
            
        gc = [data[smp]["gc"] for smp in data.keys()]
        gclist = [val for subl in gc for val in subl]
        gc_range=[min(gclist),max(gclist)]
        
        ptn = [data[smp]["ptn"].keys() for smp in data.keys()]
        ptnlist = set([val for subl in ptn for val in subl])
        
        ptn_cnt_range=[0,0]
        
        for pt in ptnlist:
            pt_instances_cnt = [data[smp]["ptn"][pt].keys() for smp in data.keys()]
            instances_cnt = set([val for subl in pt_instances_cnt for val in subl])
            if max(instances_cnt) > ptn_cnt_range[1]:
                ptn_cnt_range[1] = max(instances_cnt) 
        
        ptn_result_order = [p for p in ptnlist]
        
        header = ["Name","TotalReadCount","MappedReadCount", "DuplicateReadCount"]
        header.extend(["GC%d_%d"%(g*2,g*2+2) for g in range(gc_range[0],gc_range[1]+1)])
        
        for pt in ptn_result_order:
            header.extend(["%s_%d"%(pt,g) for g in range(ptn_cnt_range[0],ptn_cnt_range[1]+1)])
        
        outarray=[]
        outarray.append(header)
        
        for sample in data:
            result= data[sample]
            outrow = [
                      data[sample]['name'],
                      data[sample]['tol'],
                      data[sample]['map'],
                      data[sample]['dup']
                      ]
            for gc in range(gc_range[0],gc_range[1]+1):
                if data[sample]['gc'].has_key(gc):
                    outrow.append(data[sample]['gc'][gc])
                else:
                    outrow.append(self.nullvalsymbol)
            for ptn in ptn_result_order:
                for ptn_cnt in range(ptn_cnt_range[0],ptn_cnt_range[1]+1):
                    if data[sample]['ptn'][ptn].has_key(ptn_cnt):
                        outrow.append(data[sample]['ptn'][ptn][ptn_cnt])
                    else:
                        outrow.append(self.nullvalsymbol)
            outarray.append(outrow)
            
        transposed_outarray = zip(*outarray)
        util.cdm(transposed_outarray, "%s/sum_table.pickle"%outdir)
        ofh = file("%s/sum_table.csv"%outdir,'wb')
        for row in transposed_outarray:
            ofh.write(','.join([str(s) for s in row])+"\n")
        ofh.close()
        
    def _countpattern(self, ptn, dna):
        
        fw = len(re.findall(ptn,dna))
        rv = len(re.findall(util.dna_complement(ptn, reverse=True),dna))
        return fw if fw > rv else rv
            
    def _tel_perm_letter(self):
        return [
        "GTATGG",
        "GTGTAG",
        "AGTGGT",
        "GTGGAT",
        "GATGGT",
        "TGTGAG",
        "GTTGGA",
        "AGGTTG",
        "GGGATT",
        "TTAGGG"
        ]
        
        
        