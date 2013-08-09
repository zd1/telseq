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
    
__all__ = ['SeqPattern']

class SeqPattern():
        
    telomere_pattern = "TTAGGG"
    gc_interval=[0.48,0.52]
    gc_binsize=0.06
    cutoff = 7
    readlength = 100
    max_count = readlength/6
    lengthatgc =  1097938400
    
    def __init__(self, outdir, name, bamfile=None, experimental=False):
        
        self.patterns = [self.telomere_pattern]
        
        if experimental:
#             self.patterns.extend(self._tel_perm_letter()) 
            self.gc_interval=[0.40,0.60]
            self.gc_binsize=0.02
        else:
            self.gc_interval=[0.48,0.52]
            self.gc_binsize=0.06
        
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
        if not os.path.exists(os.path.join(outdir,'tmp')):
             os.makedirs(os.path.join(outdir,'tmp'))
        self.name=name
        
    def bamparser(self):
        
        # only extract flag, chromsome name, position, and the actual read
        samfields = "2,3,4,10"        
        p1 = subprocess.Popen(["samtools","view", self.bamfile], stdout=subprocess.PIPE )
        p2= subprocess.Popen(["cut","-f%s"%samfields],stdin=p1.stdout, stdout=subprocess.PIPE )
        
        counter = 0
        pct = 0.1
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
                if gcfraction >= self.gc_interval[0] and gcfraction <= self.gc_interval[1]:
                    gcbin=int(gcfraction*1.0/self.gc_binsize) # count read by its GC to every gc_binsize % GC bin
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
                
                if counter%100000==0:
                    pct =  pct+0.1 if pct+0.1<1 else 0.1
                    util.update_progress(pct, readcount=counter, showpct=False)
                counter += 1
                
        util.cdm(self.bamproperty,os.path.join(self.outdir,'tmp', "%s.bamscan.pickle"%self.name ))
#         util.cdm(self.bamproperty,"%s/%s.bamscan.pickle"%(self.outdir,self.name))
        util.update_progress(1, readcount=counter, showpct=False)
        
    @classmethod    
    def integrate(cls, outdir, ids, experimental=False, bams=None, minimum=True):
        
        nullvalsymbol="0"
        resultfiles = [os.path.join(outdir,'tmp',"%s.bamscan.pickle"%(i)) for i in ids]
        
        data={}
        dummy_name=0
        
        for pkl in resultfiles:
            result = util.cl(pkl)
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
        if experimental:
             header.extend(["GC_%d"%(g) for g in range(gc_range[0],gc_range[1]+1)])
        else:
            header.extend(["GC_%d_%d"%(SeqPattern.gc_interval[0]*100, SeqPattern.gc_interval[1]*100)])
        
        for pt in ptn_result_order:
            if experimental:
                header.extend(["%s_%d"%(pt,g) for g in range(ptn_cnt_range[0],ptn_cnt_range[1]+1)])
            else:
                header.extend(["TTAGGG_%s+"%SeqPattern.cutoff])
        
        header.append("accumulatedCount")
        header.append("lengthEstimate")
        
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
                    outrow.append(nullvalsymbol)
            for ptn in ptn_result_order:
                if experimental:
                    for ptn_cnt in range(ptn_cnt_range[0],ptn_cnt_range[1]+1):
                        if data[sample]['ptn'][ptn].has_key(ptn_cnt):
                            outrow.append(data[sample]['ptn'][ptn][ptn_cnt])
                        else:
                            outrow.append(nullvalsymbol)
                    acc = 0
                    for ptn_cnt in range(SeqPattern.cutoff, SeqPattern.max_count):
                        if data[sample]['ptn'][ptn].has_key(ptn_cnt):
                            acc += data[sample]['ptn'][ptn][ptn_cnt]
                    outrow.append(acc)
                    tl = (acc*1.0/data[sample]['gc'][gc_range[0]])*SeqPattern.lengthatgc*1.0/46000
                    outrow.append(tl)
                else:
                    acc = 0
                    for ptn_cnt in range(SeqPattern.cutoff, SeqPattern.max_count):
                        if data[sample]['ptn'][ptn].has_key(ptn_cnt):
                            acc += data[sample]['ptn'][ptn][ptn_cnt]
                    outrow.append(acc)
                    tl = (acc*1.0/data[sample]['gc'][gc_range[0]])*SeqPattern.lengthatgc*1.0/46000
                    outrow.append(tl)
            outarray.append(outrow)
        
       
        tab=""
        if len(ids) == 1:
            tab="_%s"%ids[0]
        
        util.cdm(outarray, os.path.join(outdir,"tmp","sum_table%s.pickle"%tab))
        if minimum:
            transposed_outarray = zip(*outarray)
            transposed_outarray = [transposed_outarray[0],transposed_outarray[-1]]
            outarray=zip(*transposed_outarray)
        
        ofh = file(os.path.join(outdir,"sum_table%s.csv"%tab),'wb')
        for row in outarray:
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
        
        
        