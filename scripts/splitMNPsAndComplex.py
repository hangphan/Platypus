"""
This script is used to split VCF records describing multi-nucleotide polymorphisms (MNPs)
and complex (i.e. multi-SNP) records into individual base substitutions spread across
multiple VCF records. This is done in order to compare Platypus calls with e.g. AXIOM
chip calls, as these complex variants are generally treated as multiple SNPs in the chip
calls.
"""

import sys
import gzip
import os
import re

###################################################################################################

def splitVariant(chrom, pos, theId, ref, alt, qual, filters, info, theRest):
    """
    Take information from a single record of VCF, representing a complex, multi-SNP
    variant, and return a list of strings representing the same record split into
    individual SNPs.
    """
    for index,(refBase,altBase) in enumerate(zip(ref,alt)):
        if refBase != altBase:
            yield "\t".join([chrom, str(pos+index), theId, refBase, altBase, qual, filters, "%s;%s" %(info,"FromComplex"), theRest])

###################################################################################################

def splitMAVariant(chrom, pos, theId, ref, alts, qual, filters, info, theRest):
    """
    Take information from a single record of VCF, representing a complex, multi-SNP
    variant, and return a list of strings representing the same record split into
    individual SNPs.
    """
    
    splitVars = {} # Store pos, ref, alt
    lenVar = len(ref)
    for index in range(lenVar):
        refBase = ref[index]
        altBases = []
        for alt in alts:
            altBases.append(alt[index])
        if len(set(altBases)) ==1 and refBase in set(altBases):
            continue
        splitVars[(pos + index, refBase)] = altBases
                    
    cols = theRest.split()
    for (thisPos, refBase) in sorted(splitVars.keys()):
        allUniqAltBases = [refBase] + sorted(list(set(splitVars[(thisPos, refBase)]).difference(set(refBase))))
        allAltBases = splitVars[(thisPos, refBase)]
        splitVars[(thisPos, refBase)] = allUniqAltBases[1:]
        newRest = [cols[0]]
        for sample in cols[1:]:
            newField = sample.split(":")
            gt1 = newField[0][0]
            gt2 = newField[0][2]
            altBases = ['', '']
            if gt1 != '.':
                gt1 = int(newField[0][0]) -1
                altBases[0] = allAltBases[gt1]
                gt1 = allUniqAltBases.index(altBases[0]) 
            if gt2 != '.':
                gt2 = int(newField[0][2]) -1
                altBases[1] = allAltBases[gt2]
                gt2 = allUniqAltBases.index(altBases[1]) 
            newField[0] = '/'.join(map(str,[gt1, gt2]))
            newRest.append( ":".join(newField))
        yield "\t".join([chrom, str(thisPos), theId, refBase, ','.join(splitVars[(thisPos, refBase)]), qual, filters, "%s;%s" %(info,"FromComplex"), '\t'.join(newRest)])



###################################################################################################

for line in sys.stdin:

    if line.startswith("#"):
        print line.strip()
        continue
    else:
        col = line.strip().split("\t")
        chrom = col[0]
        pos = int (col[1])
        theId = col[2]
        ref = col[3]
        alts = col[4].split(",")
        qual = col[5]
        filters = col[6]
        info = col[7]
        theRest = "\t".join(col[8:])

        # Multi-allelic sites and length-changing variants are not further processed,
        # and go straight to output. Also, don't split complex variants/MNPs of length < 4
        if len(ref) != len(alts[0]) or len(ref) < 2 or (len(alts) > 1 and len(alts[1]) != len(ref)) or (len(alts) > 2 and len(alts[2]) != len(ref)) or len(alts) > 3:
            print line.strip()
        # Special treatment for multi-allelic sites
        elif len(alts) > 1:
            for newLine in splitMAVariant(chrom, pos, theId, ref, alts, qual, filters, info, theRest):
                print newLine
        # Otherwise, split the variant into SNPs
        else:
            for newLine in splitVariant(chrom, pos, theId, ref, alts[0], qual, filters, info, theRest):
                print newLine
