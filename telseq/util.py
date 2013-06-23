'''
-----------------------------------------------
Copyright 2013 Wellcome Trust Sanger Institute
Written by Zhihao Ding (zd1@sanger.ac.uk)
Released under the GPL
-----------------------------------------------

Created on 26 May 2013

'''

import re
import cPickle

def cl(f): 
    return cPickle.load(open(f,'rb'))
def cdm(o,f): 
    cPickle.dump(o, open(f, 'wb'), -1)

def dna_complement(dna, reverse=False):
    
    rvs_dna=[]
    nucleotide_mapping = {
        'A':'T',
        'G':'C',
        'T':'A',
        'C':'G'
                          }
    for letter in dna:
        if nucleotide_mapping.has_key(letter.upper()):
            rvs_dna.append(nucleotide_mapping[letter.upper()])
        else:
            rvs_dna.append("N")
            
    if reverse:
        return "".join(rvs_dna[::-1])
    else:
        return "".join(rvs_dna)

def gc_fraction(dna):
    return len(re.findall('[GCgc]', dna))*1.0/len(dna)
    
def flagresolve(flagnumber):
    
    flagnumber = int(flagnumber)
    bitflag = {
            0x1:'template having multiple segments in sequencing',
            0x2:"each segment properly aligned according to the aligner", 
            0x4: "segment unmapped",
            0x8: "next segment in the template unmapped",
            0x10: "SEQ being reverse complemented",
            0x20: "SEQ of the next segment in the template being reversed", 
            0x40: "the first segment in the template",
            0x80: "the last segment in the template",
            0x100: "secondary alignment",
            0x200: "not passing quality controls", 
            0x400: "PCR or optical duplicate"
    }
    
    setflags = {}
    for key in bitflag.keys():
        if key & flagnumber:
            setflags[key] = bitflag[key]
    
    return setflags
            

