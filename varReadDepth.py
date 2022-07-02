#!/usr/bin/env python
'''
	Author: Ti-Cheng Chang
	Purporse: Extract the cell barcode and UMI index , and counts for single cell data
'''



import pysam, re, argparse, math, logging
import collections
import pandas as pd
from collections import OrderedDict



logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.INFO,
    format="[%(levelname)s][%(asctime)s]%(message)s",
    handlers=[
        logging.FileHandler('run_varReadDepth.log', mode='w'),
        logging.StreamHandler()
    ]
)



def _var_read_depth(cds, bam, output, id, tag1='CB', tag2='UB', minimum_base_quality = 15, minimum_mapq = 255, minimum_read_quality = 15, verbose=5, ALT_n_LB=None, write_output=True, dedup=False, print_reads=False, print_base=None, print_pos=None, stepper='all'):
    skip_message = 'SKIPPED: {} {} {} {} {}'
    pass_message = 'PASSED: {} {} {} {}'
    out=[]

    var=pd.read_csv(cds, sep='\t', index_col=False)
    var.columns = ['chr','position','ref','alt']
    print(var)
    print('variation table size', 'nrow', var.shape[0], 'ncol', var.shape[1])
    var.drop_duplicates(inplace=True)
    print('variation table size after filter duplicates', 'nrow', var.shape[0], 'ncol', var.shape[1]) 
    var.dropna(axis=0, how='any',inplace=True)
    print('variation table size after filter NA', 'nrow', var.shape[0], 'ncol', var.shape[1]) 

    header=['ID', 'chr', 'position', 'total_n', 'QC_n', 'ref', 'alt', tag1 , tag2, 'refcount', 'altcount', 
            'good_readN', 'good_A', 'good_T', 'good_C', 'good_G', 
            'ins_readN', 'del_readN', 
            'lowmapq_readN', 'lowmapq_A', 'lowmapq_T', 'lowmapq_C', 'lowmapq_G',
            'nomapq_readN', 'nomapq_A', 'nomapq_T', 'nomapq_C', 'nomapq_G',
            'lowreadq_readN', 'lowreadq_A', 'lowreadq_T', 'lowreadq_C', 'lowreadq_G',
            'lowbaseq_readN', 'lowbaseq_A', 'lowbaseq_T', 'lowbaseq_C', 'lowbaseq_G']   
    fh=open(output, 'w')
    #print '\t'.join(header)
    fh.write('\t'.join(header) + '\n')

    target_base=list()
    if print_reads:
        output2=output+'.reads.txt'
        fh2=open(output2, 'w')
        if print_base is not None:
            target_base=print_base.split(',') 
        

    for idx, row in var.iterrows():
        refbase=row['ref']
        altbase=row['alt']
        samfile = pysam.AlignmentFile(bam, mode="rb") #, ignore_truncation=False) 

        trimchr=False
        if not any('chr' in s for s in samfile.references):trimchr=True
        if trimchr:
            tarchr=str(row['chr']).lstrip('chr') 
        else:
            tarchr=str(row['chr'])
            if 'chr' not in tarchr: tarchr='chr'+tarchr
    

        #print(tarchr, row['position'], refbase, altbase)
        for column in samfile.pileup(tarchr, int(row['position'])-1, int(row['position']), max_depth=1000000000, stepper=stepper):
            pos0based=column.pos
            pos1based=column.pos+1
            if int(row['position']) == int(pos1based):
                 barcode         ={}
                 lowmapq_barcode ={}
                 nomapq_barcode  ={}
                 lowreadq_barcode={}
                 lowbaseq_barcode={}

                 good_count    ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 ins_count     ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 del_count     ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 lowmapq_count ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 nomapq_count  ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 lowreadq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 lowbaseq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 
                 ### process each read in each exon to get the count of bases
                 good_n=0
                 del_n =0
                 ins_n =0
                 for pileupread in column.pileups:
                     read     = pileupread.alignment
                     size     = pileupread.indel
                     qpos     = pileupread.query_position
                     qname    = pileupread.alignment.query_name
                    
                     Tag1='no'+tag1
                     if read.has_tag(tag1):
                         Tag1=read.get_tag(tag1) ### barcode: CB, 

                     Tag2='no'+tag2
                     if read.has_tag(tag2):
                         Tag2=read.get_tag(tag2) ### barcode: CB,

                     if dedup and read.is_duplicate: continue


                     if refbase=='-' or len(altbase) > len(refbase):
                         ### INSERTION
                         ins_len=len(altbase)
                         basecount=0
                         prev_is=None
                         prev_pos=None
                         inseq=''
                         try:
                          basecomposition=read.get_aligned_pairs(with_seq=True, matches_only=True)
                          for b in basecomposition:
                            if pos1based == b[1]+1:
                                basecount+=1
                            if pos1based <= b[1]+1 <= pos1based+ins_len+1:
                                if prev_is is None:
                                    prev_is=b[0]
                                    prev_pos=b[1]+1
                                else:
                                    if b[0] - prev_is > 1: ### with insertion
                                        #print 'INSpos:', prev_pos, b[1]+1, prev_is+1, b[0]-1
                                        for i in range((prev_is+1), (b[0])):
                                            inseq+=read.query_sequence[i]
                                    prev_is=b[0]
                                    prev_pos=b[1]+1
                         except:
                             bascount=-1

                         if basecount > 0 and Tag1 not in barcode:
                             barcode[Tag1]={}
                             barcode[Tag1][Tag2]={'INS':0, 'noINS':0}
                             barcode[Tag1][Tag2]={'DEL':0, 'noDEL':0}
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowmapq_barcode[Tag1]={}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1]={}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1]={}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1]={}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

                         if basecount > 0 and Tag2 not in barcode[Tag1]:
                             barcode[Tag1][Tag2]={'INS':0, 'noINS':0}
                             barcode[Tag1][Tag2]={'DEL':0, 'noDEL':0}
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

    
                         #print '[INS with ins]', row['chr'], pos1based, pos1based+ins_len, refbase, qname, inseq
 
                         if len(inseq) > 0 and str(len(inseq))+'I' in read.cigarstring: 
                             #print '[INS with ins]', row['chr'], pos1based, pos1based+ins_len, refbase, qname, inseq
                             barcode[Tag1][Tag2]['INS']+=1   
                         elif basecount > 0:
                             barcode[Tag1][Tag2]['noINS']+=1

                         if basecount>0:
                             good_n+=1

                         if print_reads:
                             fh2.write('\t'.join(list(map(str,[qname, pos1based,'INS'])))+'\n')

                     elif altbase=='-' or len(refbase) > len(altbase):
                         ### DELETION
                         del_len=len(refbase)
                         #print '[DEL]', row['chr'], pos1based, pos1based+del_len, refbase
                         basecount=0
                         try:
                          basecomposition=read.get_aligned_pairs(with_seq=True, matches_only=True)
                          for b in basecomposition:
                            if pos1based-1 <= b[1]+1 <= pos1based+del_len:
                                print("TEST")
                                print(b)
                                basecount+=1
                                print('[DEL]', row['chr'], pos1based, pos1based+del_len, refbase, qname, b[1]+1, b[2])
                         except:
                             basecount=-1                     

                         if basecount > 0 and Tag1 not in barcode:
                             barcode[Tag1]={}
                             barcode[Tag1][Tag2]={'INS':0, 'noINS':0}
                             barcode[Tag1][Tag2]={'DEL':0, 'noDEL':0}
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowmapq_barcode[Tag1]={}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1]={}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1]={}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1]={}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

                         if basecount > 0 and Tag2 not in barcode[Tag1]:
                             barcode[Tag1][Tag2]={'INS':0, 'noINS':0}
                             barcode[Tag1][Tag2]={'DEL':0, 'noDEL':0}
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}


                         if basecount==2 and (str(del_len)+'D' in read.cigarstring): ## or str(del_len)+'N' in read.cigarstring):
                             #print '[DEL]', row['chr'], pos1based, pos1based+del_len, refbase, qname
                             barcode[Tag1][Tag2]['DEL']+=1
                         elif basecount > 0:
                             barcode[Tag1][Tag2]['noDEL']+=1
                         
                         if basecount > 0:
                             good_n+=1       

                         if print_reads:
                             fh2.write('\t'.join(list(map(str,[qname, pos1based,'DEL'])))+'\n')

 
                     else:
                         ### SNV
                         #print read
                         if pileupread.is_del: #or size < 0:
                             readbase=''
                         else:
                             readbase = read.query_sequence[qpos]

                         #if read.is_duplicate:
                         #    print '[dup]', row['chr'], pos1based, qname, size, qpos, read.cigarstring, readbase
                         # print("\nregion:%s:%s-%s, coverage at base %s:%s = %s" % (row['chr'], row['position']-1, row['position'], pos1based, refbase, column.n))
                         if readbase not in ['A', 'T', 'C', 'G']:
                             continue

                         #print row['chr'], pos1based, qname, size, qpos, read.cigarstring, readbase
                     
                         if Tag1 not in barcode:
                             barcode[Tag1]={}
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             barcode[Tag1][Tag2]={'INS':0, 'noINS':0}
                             barcode[Tag1][Tag2]={'DEL':0, 'noDEL':0}
                             lowmapq_barcode[Tag1]={}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1]={}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1]={}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1]={}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

                         if Tag2 not in barcode[Tag1]:
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             barcode[Tag1][Tag2]={'INS':0, 'noINS':0}
                             barcode[Tag1][Tag2]={'DEL':0, 'noDEL':0}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                         #print pos1based, qname, size, qpos, read.cigarstring, readbase, Tag
                     

                         # skip reads with indels
                         if pileupread.is_del:# or pileupread.indel < 0:
                             print(skip_message.format(pos1based, qname, 'is del', readbase, pileupread.indel))
                             del_count[refbase]+=1
                             del_n+=1
                             continue

                         ## skip reads with mapq below threshold
                         if minimum_mapq is not None and pileupread.alignment.mapping_quality < minimum_mapq:
                         #    print skip_message.format(pos1based, qname, 'low mapq', readbase, pileupread.alignment.mapping_quality)
                             lowmapq_count[readbase]+=1
                             lowmapq_barcode[Tag1][Tag2][readbase]+=1
                             continue
                         ## skip reads with no mapq specified
                         ## Star uses a similar scoring scheme to tophat except that the value for uniquely mapped reads is 255 instead of 50.
                         ## The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10*log10(1-1/[number of loci the read maps to])) for multi-mapping reads. This scheme is same as the one used by Tophat.
                         elif read.mapping_quality is None:
                         #    print pass_message.format(pos1based, qname, 'no mapq', '255')
                             nomapq_count[readbase]+=1
                             nomapq_barcode[Tag1][Tag2][readbase]+=1
                             continue
                         ## skip mean qscore of the read below threshold
                         elif mean(read.query_qualities) < minimum_read_quality:
                         #    print pass_message.format(pos1based, qname, 'low read quality', mean(read.query_qualities))
                             lowreadq_count[readbase]+=1
                             lowmapq_barcode[Tag1][Tag2][readbase]+=1
                             continue
                         else:
                             # check for insertion
                             if pileupread.is_refskip: #or pileupread.indel > 0:
                                 if read.query_qualities[qpos+1] >= minimum_base_quality:
                                      insbase=read.query_sequence[qpos+1]
                                      insbase=complement(insbase) if strand =='-' else insbase
                                      #print skip_message.format(pos1based, qname, 'is ins', readbase)
                                      ins_count[insbase]+=1
                                      ins_n+=1
                             # skip reads with a base quality below threshold
                             if read.query_qualities[qpos] < minimum_base_quality:
                                  #print skip_message.format(position, qname, 'low base quality', read.query_qualities[pileupread.query_position])
                                  lowbaseq_count[readbase]+=1
                                  lowbaseq_barcode[Tag1][Tag2][readbase]+=1
                                  continue
 
                         ### print details of each passed read
                         #print pass_message.format(position, qname, read.query_sequence[qpos])
                         good_n+=1
                         good_count[readbase]+=1
                         print(Tag1, Tag2, readbase)
                         print(barcode)
                         barcode[Tag1][Tag2][readbase]+=1
                         if print_reads:
                             if readbase in target_base:
                                 fh2.write('\t'.join(list(map(str,[qname, pos1based, readbase])))+'\n')

                 #print '[good]', good_n, refbase

                 good_readN    =readNSum(good_count)
                 ins_readN     =readNSum(ins_count)
                 del_readN     =readNSum(del_count)
                 lowmapq_readN =readNSum(lowmapq_count)
                 nomapq_readN  =readNSum(nomapq_count)
                 lowreadq_readN=readNSum(lowreadq_count)
                 lowbaseq_readN=readNSum(lowbaseq_count)
 
 
                 #if good_n > 0:
                 if refbase not in ['A','T','C','G','N']:
                     for b in barcode:
                             for u in barcode[b]:
                                 out=map(str, [id, row['chr'], pos1based, column.n, good_n, refbase, altbase, b , u, 'NA', 'NA',
                                               good_readN, barcode[b][u]['A'], barcode[b][u]['T'], barcode[b][u]['C'], barcode[b][u]['G'],
                                               ins_readN, del_readN,
                                               lowmapq_readN, lowmapq_barcode[b][u]['A'], lowmapq_barcode[b][u]['T'], lowmapq_barcode[b][u]['C'], lowmapq_barcode[b][u]['G'],
                                               nomapq_readN, nomapq_barcode[b][u]['A'], nomapq_barcode[b][u]['T'], nomapq_barcode[b][u]['C'], nomapq_barcode[b][u]['G'],
                                               lowreadq_readN, lowreadq_barcode[b][u]['A'], lowreadq_barcode[b][u]['T'], lowreadq_barcode[b][u]['C'], lowreadq_barcode[b][u]['G'],
                                               lowbaseq_readN, lowbaseq_barcode[b][u]['A'], lowbaseq_barcode[b][u]['T'], lowbaseq_barcode[b][u]['C'], lowbaseq_barcode[b][u]['G'] ])
                                 #print '\t'.join(out)
                                 fh.write('\t'.join(out) + '\n')


                 
                 elif len(refbase)==1 and refbase != '-' and len(altbase)==1:
                         #print 'SNV'
                         for b in barcode:
                             for u in barcode[b]:
                                 out=map(str, [id, row['chr'], pos1based, column.n, good_n, refbase, altbase, b , u, barcode[b][u][refbase], barcode[b][u][altbase], 
                                               good_readN, barcode[b][u]['A'], barcode[b][u]['T'], barcode[b][u]['C'], barcode[b][u]['G'], 
                                               ins_readN, del_readN, 
                                               lowmapq_readN, lowmapq_barcode[b][u]['A'], lowmapq_barcode[b][u]['T'], lowmapq_barcode[b][u]['C'], lowmapq_barcode[b][u]['G'],
                                               nomapq_readN, nomapq_barcode[b][u]['A'], nomapq_barcode[b][u]['T'], nomapq_barcode[b][u]['C'], nomapq_barcode[b][u]['G'],
                                               lowreadq_readN, lowreadq_barcode[b][u]['A'], lowreadq_barcode[b][u]['T'], lowreadq_barcode[b][u]['C'], lowreadq_barcode[b][u]['G'],
                                               lowbaseq_readN, lowbaseq_barcode[b][u]['A'], lowbaseq_barcode[b][u]['T'], lowbaseq_barcode[b][u]['C'], lowbaseq_barcode[b][u]['G'] ])
                                      
                                 #print '\t'.join(out)
                                 fh.write('\t'.join(out) + '\n')
                 
                 elif altbase == '-':
                         #print 'DEL'
                         for b in barcode:
                             for u in barcode[b]:
                                 out=map(str, [id, row['chr'], pos1based, column.n, good_n, refbase, altbase, b , u, barcode[b][u]['noDEL'], barcode[b][u]['DEL'], 
                                               good_readN, 'NA', 'NA', 'NA', 'NA',
                                               ins_readN, del_readN, 
                                               lowmapq_readN, 'NA', 'NA', 'NA', 'NA',
                                               nomapq_readN, 'NA', 'NA', 'NA', 'NA',
                                               lowreadq_readN, 'NA', 'NA', 'NA', 'NA',
                                               lowbaseq_readN, 'NA', 'NA', 'NA', 'NA'] )
                                 #print '\t'.join(out)
                                 fh.write('\t'.join(out) + '\n')
                 elif refbase == '-':
                         #print 'INS'
                         for b in barcode:
                             for u in barcode[b]:
                                 out=map(str, [id, row['chr'], pos1based, column.n, good_n, refbase, altbase, b , u, barcode[b][u]['noINS'], barcode[b][u]['INS'], 
                                               good_readN, 'NA', 'NA', 'NA', 'NA', 
                                               ins_readN, del_readN, 
                                               lowmapq_readN, 'NA', 'NA', 'NA', 'NA',
                                               nomapq_readN, 'NA', 'NA', 'NA', 'NA',
                                               lowreadq_readN, 'NA', 'NA', 'NA', 'NA',
                                               lowbaseq_readN, 'NA', 'NA', 'NA', 'NA'] )
                                 #print '\t'.join(out)
                                 fh.write('\t'.join(out) + '\n')
                 else:
                         #print 'Others'
                         out=map(str, [id, row['chr'], pos1based, column.n, good_n, refbase, altbase, 'NA' , 'NA', 'NA', 'NA', 
                                       good_readN, 'NA', 'NA', 'NA', 'NA',
                                       ins_readN, del_readN, 
                                       lowmapq_readN, 'NA', 'NA', 'NA', 'NA',
                                       nomapq_readN, 'NA', 'NA', 'NA', 'NA',
                                       lowreadq_readN, 'NA', 'NA', 'NA', 'NA',
                                       lowbaseq_readN, 'NA', 'NA', 'NA', 'NA', ])
                         #print '\t'.join(out) 
                         fh.write('\t'.join(out) + '\n')
    fh.close()
    fh2.close()


def var_read_depth(cds, bam, output, id, reffn, tag1='CB', tag2='UB', \
                   minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, max_d=999999999, dedup=True, \
                   verbose=5, ALT_n_LB=None, write_output=True, \
                   print_reads=False, print_base=None, print_pos=None, stepper='no_filter',skip_indel=False, \
                   barcode_reads=None, drop_na=False):
    '''
        NOTE: need to add dedup function
    '''
    #max_d=200 ### test
    skip_message = 'SKIPPED: {} {} {} {} {} {}'
    pass_message = 'PASSED: {} {} {} {} {}'
    out=[]
  
    var=pd.read_csv(cds, sep='\t', index_col=False)
    var.columns = ['chr','position','ref','alt']
    print(var)
    print('variation table size', 'nrow', var.shape[0], 'ncol', var.shape[1])
    var.drop_duplicates(inplace=True)
    print('variation table size after filter duplicates', 'nrow', var.shape[0], 'ncol', var.shape[1])
    if drop_na:
        var.dropna(axis=0, how='any',inplace=True)
        print('variation table size after filter NA', 'nrow', var.shape[0], 'ncol', var.shape[1])

    ref=pysam.FastaFile(reffn)

    header=['ID', 'chr', 'position', 'total_n', 'QC_n_no_overlap', 'ref', 'alt', tag1 , tag2, 'refcount', 'altcount',
            'good_readN', 'good_A', 'good_T', 'good_C', 'good_G',
            #'good_Af', 'good_Tf', 'good_Cf', 'good_Gf',
            #'good_Ar', 'good_Tr', 'good_Cr', 'good_Gr',
            'ins_readN', 'del_readN',
            'lowmapq_readN', 'lowmapq_A', 'lowmapq_T', 'lowmapq_C', 'lowmapq_G',
            'nomapq_readN', 'nomapq_A', 'nomapq_T', 'nomapq_C', 'lowmapq_G',
            'lowreadq_readN', 'lowreadq_A', 'lowreadq_T', 'lowreadq_C', 'lowreadq_G',
            'lowbaseq_readN', 'lowbaseq_A', 'lowbaseq_T', 'lowbaseq_C', 'lowbaseq_G']
    if write_output:
        fh=open(output, 'w')
        fh.write('\t'.join(header) + '\n')
    #else:
    #    print('\t'.join(header))

    target_base=list()
    if print_reads:
        output2=output+'.reads.txt'
        fh2=open(output2, 'w')
        fh2.write('\t'.join(list(map(str,['read', 'pos', 'base'])))+'\n')
        if print_base is not None:
            target_base=print_base.split(',')


    samfile = pysam.AlignmentFile(bam, mode="rb") #, ignore_truncation=False) 

    trimchr=False
    if not any('chr' in s for s in samfile.references):trimchr=True

    #logging.debug('[BAM]' + str(bam))
   
    for idx, row in var.iterrows():
        chr=row['chr']
        pos=int(row['position'])
        refbase=row['ref']
        altbase=row['alt']
        tarchr=str(chr).lstrip('chr') if trimchr else chr
        
        print(tarchr, pos, refbase, altbase)
 
        for column in samfile.pileup(tarchr, pos-1, pos, max_depth=max_d, stepper=stepper, truncate=True, ignore_overlaps=False, min_base_quality=0):
            pos0based=column.pos
            pos1based=column.pos+1
            if int(pos) == int(pos1based):
                 #logging.debug('[Var depth] pos found:' + str(pos) + '; depth:' + str(column.n))
                 barcode         ={}
                 lowmapq_barcode ={}
                 nomapq_barcode  ={}
                 lowreadq_barcode={}
                 lowbaseq_barcode={}

                 good_count    ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 good_count_reads={'A':list(), 'T':list(), 'C':list(), 'G':list(), 'N':list()}
                 good_fcount   ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 good_rcount   ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 ins_count     ={'Y':0, 'N':0}
                 ins_fcount    ={'Y':0, 'N':0}
                 ins_rcount    ={'Y':0, 'N':0}
                 del_count     ={'Y':0, 'N':0}
                 del_fcount    ={'Y':0, 'N':0}
                 del_rcount    ={'Y':0, 'N':0}
                 mm_count     ={'Y':0, 'N':0}
                 mm_fcount    ={'Y':0, 'N':0}
                 mm_rcount    ={'Y':0, 'N':0}
                 lowmapq_count ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 nomapq_count  ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 lowreadq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                 lowbaseq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

                 ### process each read in each exon to get the count of bases
                 good_n=0
                 del_n =0
                 ins_n =0
                 #for tst in column.pileups:
                 #    print(tst.alignment.query_name)

                 for pileupread in column.pileups:
                     read     = pileupread.alignment
                     size     = pileupread.indel
                     qpos     = pileupread.query_position
                     qname    = pileupread.alignment.query_name
                      

                     if dedup and read.is_duplicate:
                         continue
                     #if read.is_qcfail or read.is_secondary or read.is_supplementary or read.is_unmapped:
                     #    continue

                     Tag1='no'+tag1
                     if read.has_tag(tag1):
                         Tag1=read.get_tag(tag1) ### barcode: CB, 

                     Tag2='no'+tag2
                     if read.has_tag(tag2):
                         Tag2=read.get_tag(tag2) ### barcode: CB,


                     if refbase=='-' or (len(altbase) > 1 and len(refbase)==1):
                       if skip_indel: break
                       ### INSERTION
                       if not read.has_tag('MD'):
                           logging.debug('[INS] No MD tag')
                           break
                       ins_len=len(altbase)-len(refbase)
                       try:
                          basecomposition=read.get_aligned_pairs(with_seq=True, matches_only=True)
                          #print('INS', tarchr, pos, refbase, altbase)
                       except:
                          break
                       basecount=0
                       prev_is=None
                       prev_pos=None
                       inseq=''
                       read_ins_start=None
                       read_ins_end=None
                       for b in basecomposition:
                            if pos1based == b[1]+1:
                                read_ins_start=b[0]
                            if pos1based+1 == b[1]+1:
                                read_ins_end=b[0]
                                #logging.info('[INS]', qname, read_ins_start, read_ins_end, read.seq[read_ins_start:read_ins_end])
                                if read.seq[read_ins_start:read_ins_end] == altbase:
                                      ins_count['Y']+=1
                                      if read.is_reverse:
                                          ins_rcount['Y']+=1
                                      else:
                                          ins_fcount['Y']+=1
                                      if print_reads:
                                          fh2.write('\t'.join(list(map(str,[qname, pos1based, read.seq[read_ins_start:read_ins_end]])))+'\n')
                                else:
                                      ins_count['N']+=1
                                      if read.is_reverse:
                                          ins_rcount['N']+=1
                                      else:
                                          ins_fcount['N']+=1
                                      if print_reads:
                                          fh2.write('\t'.join(list(map(str,[qname, pos1based, refbase])))+'\n')
                                break
                     elif altbase=='-' or (len(refbase) > 1 and len(altbase)==1):
                         if skip_indel: break
                         ### DELETION
                         del_len=len(refbase)-len(altbase)
                         #print '[DEL]', row['chr'], pos1based, pos1based+del_len, refbase
                         if not read.has_tag('MD'):
                           logging.debug('[DEL] No MD tag')
                           break
                         try:
                           basecomposition=read.get_aligned_pairs(with_seq=True, matches_only=True)
                         except: break
                         basecount=0
                         read_del_start=None
                         ref_del_start=None
                         ref_del_end=None
                         for b in basecomposition:
                            if pos1based == b[1]+1:
                                ref_del_start=b[1]
                                read_del_start=b[0]
                            if read_del_start is not None and b[0] == read_del_start+1:
                                ref_del_end=b[1]
                                del_ref=ref.fetch(tarchr, ref_del_start, ref_del_end)
                                #logging.info('[DEL]', qname, ref_del_start, ref_del_end, del_ref)
                                if del_ref==refbase:
                                    del_count['Y']+=1
                                    if read.is_reverse:
                                        del_rcount['Y']+=1
                                    else:
                                        del_fcount['Y']+=1
                                    if print_reads:
                                        fh2.write('\t'.join(list(map(str,[qname, pos1based, altbase])))+'\n')
                                else:
                                    del_count['N']+=1
                                    if read.is_reverse:
                                       del_rcount['N']+=1
                                    else:
                                       del_fcount['N']+=1
                                    if print_reads:
                                        fh2.write('\t'.join(list(map(str,[qname, pos1based, refbase])))+'\n')
                                
                                break
                     elif len(refbase) > 1 and len(altbase) > 1:
                             if skip_indel: break
                             if not read.has_tag('MD'):
                               logging.debug('[Complex] No MD tag')
                               break
                             try:
                                basecomposition=read.get_aligned_pairs(with_seq=True, matches_only=True)
                             except: break
                             read_mm_start=None
                             read_mm_end=None
                             ref_mm_start=None
                             ref_mm_end=None
                             read_mm=None
                             ref_mm=None
                             for b in basecomposition:
                                 if pos1based==b[1]+1:
                                     ref_mm_start=b[1]
                                     read_mm_start=b[0]
                                 if read_mm_start is not None and b[0]==read_mm_start+len(refbase):
                                     read_mm_end=b[0]
                                     read_mm=read.seq[read_mm_start:read_mm_end]
                                 if ref_mm_start is not None and b[1]==ref_mm_start+len(altbase):
                                     ref_mm_end=ref_mm_start+len(altbase)
                                     ref_mm=ref.fetch(tarchr, ref_mm_start, ref_mm_end)
                             if refbase==ref_mm and altbase==read_mm:
                                  mm_count['Y']+=1
                                  if read.is_reverse:
                                        mm_rcount['Y']+=1
                                  else:
                                        mm_fcount['Y']+=1
                                  if print_reads:
                                        fh2.write('\t'.join(list(map(str,[qname, pos1based, read_mm])))+'\n')
                             else:
                                  mm_count['N']+=1
                                  if read.is_reverse:
                                       mm_rcount['N']+=1
                                  else:
                                       mm_fcount['N']+=1
                                  if print_reads:
                                        fh2.write('\t'.join(list(map(str,[qname, pos1based, ref_mm])))+'\n')
                     else:
                         ### SNV
                         if pileupread.is_del: #or size < 0:
                             readbase=''
                         else:
                             readbase = read.query_sequence[qpos]
                         #if read.is_duplicate:
                         #    print '[dup]', row['chr'], pos1based, qname, size, qpos, read.cigarstring, readbase
                         # print("\nregion:%s:%s-%s, coverage at base %s:%s = %s" % (row['chr'], row['position']-1, row['position'], pos1based, refbase, column.n))
                         if readbase not in ['A', 'T', 'C', 'G']:
                             #! logging.debug(' '.join(map(str, ['[Var depth]', tarchr, pos1based, qname, 'unknown base:', readbase])))
                             continue
                             #return {'A':'.', 'T':'.', 'C':'.', 'G':'.', 'N':'.'},{'A':'.', 'T':'.', 'C':'.', 'G':'.', 'N':'.'},{'A':'.', 'T':'.', 'C':'.', 'G':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'},{'Y':'.', 'N':'.'}
                         else:
                             ori='R' if read.is_reverse else 'F'
                             #! logging.debug(' '.join(map(str,['[Var depth]',tarchr, pos1based, qname, 'base:', readbase, 'orientation:', ori])))
                         #print(row['chr'], pos1based, qname, size, qpos, read.cigarstring, readbase)


                         if Tag1 not in barcode:
                             barcode[Tag1]={}
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowmapq_barcode[Tag1]={}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1]={}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1]={}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1]={}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

                         if Tag2 not in barcode[Tag1]:
                             barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowmapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             nomapq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowreadq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                             lowbaseq_barcode[Tag1][Tag2]={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                         #print pos1based, qname, size, qpos, read.cigarstring, readbase, Tag

                         # skip reads with indels
                         #if pileupread.is_del:# or pileupread.indel < 0:
                         #    #print skip_message.format(pos1based, qname, 'is del', readbase, pileupread.indel)
                         #    del_count[refbase]+=1
                         #    del_n+=1
                         #    continue

                         ## skip reads with mapq below threshold
                         if minimum_mapq is not None and pileupread.alignment.mapping_quality < minimum_mapq:
                             #! logging.debug(' '.join(list(map(str, ['[low mapq]', tarchr, pos1based, qname, readbase, pileupread.alignment.mapping_quality]))))
                             lowmapq_count[readbase]+=1
                             lowmapq_barcode[Tag1][Tag2][readbase]+=1
                             continue
                         ## skip reads with no mapq specified
                         ## Star uses a similar scoring scheme to tophat except that the value for uniquely mapped reads is 255 instead of 50.
                         ## The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10*log10(1-1/[number of loci the read maps to])) for multi-mapping reads. This scheme is same as the one used by Tophat.
                         elif read.mapping_quality is None:
                             #! logging.debug(' '.join(list(map(str, ['[no mapq]', pos1based, qname, 'no mapq', '255']))))
                             nomapq_count[readbase]+=1
                             nomapq_barcode[Tag1][Tag2][readbase]+=1
                             continue
                         ## skip mean qscore of the read below threshold
                         elif mean(read.query_qualities) < minimum_read_quality:
                             #! logging.debug(' '.join(list(map(str, ['[low read quality]', tarchr, pos1based, qname, mean(read.query_qualities)]))))
                             lowreadq_count[readbase]+=1
                             lowmapq_barcode[Tag1][Tag2][readbase]+=1
                             continue
                         else:
                             # check for insertion
                             if pileupread.is_refskip: #or pileupread.indel > 0:
                                 if read.query_qualities[qpos+1] >= minimum_base_quality:
                                      insbase=read.query_sequence[qpos+1]
                                      insbase=complement(insbase) if strand =='-' else insbase
                                      #print skip_message.format(pos1based, qname, 'is ins', readbase)
                                      ins_count[insbase]+=1
                                      ins_n+=1
                             # skip reads with a base quality below threshold
                             if read.query_qualities[qpos] < minimum_base_quality:
                                  msg=['[low base quality]', tarchr, pos1based, qname, read.query_qualities[pileupread.query_position]]
                                  #! logging.debug(' '.join(list(map(str, msg))))
                                  lowbaseq_count[readbase]+=1
                                  lowbaseq_barcode[Tag1][Tag2][readbase]+=1
                                  continue
                         ### print details of each passed read
                         #print(pass_message.format(position, qname, read.query_sequence[qpos]))
                         good_n+=1
                         good_count[readbase]+=1
                         if qname not in good_count_reads[readbase]:
                             good_count_reads[readbase].append(qname)
                         barcode[Tag1][Tag2][readbase]+=1
                         ### count strand
                         if read.is_reverse:
                             good_rcount[readbase]+=1
                         else:
                             good_fcount[readbase]+=1

                         if print_reads:
                             if len(target_base)> 0 and readbase in target_base:
                                 fh2.write('\t'.join(list(map(str,[qname, pos1based, readbase])))+'\n')
                             else:
                                 fh2.write('\t'.join(list(map(str,[qname, pos1based, readbase])))+'\n')


                 ins_readN     =readNSum(ins_count)
                 del_readN     =readNSum(del_count)
                 lowmapq_readN =readNSum(lowmapq_count)
                 nomapq_readN  =readNSum(nomapq_count)
                 lowreadq_readN=readNSum(lowreadq_count)
                 lowbaseq_readN=readNSum(lowbaseq_count)
                 good_readN    =readNSum(good_count) + ins_readN + del_readN
                 good_count_no_overlap={'A':len(good_count_reads['A']),
                                        'T':len(good_count_reads['T']),
                                        'C':len(good_count_reads['C']),
                                        'G':len(good_count_reads['G']),
                                        'N':len(good_count_reads['N'])}
                 good_count_no_overlapN=readNSum(good_count_no_overlap) + ins_readN + del_readN
                 #print(ins_count, del_count, mm_count)
                 #! logging.debug(' '.join(map(str,["[var_read_depth]",chr,pos,good_count_no_overlap, good_fcount, good_rcount, ins_count, ins_fcount, ins_rcount, del_count, del_fcount, del_rcount, mm_count, mm_fcount, mm_rcount])))
                 #del samfile, good_count_reads
                 #gc.collect()
                 if write_output:
                        refcount=-1
                        altcount=-1
                        #if refbase in ['A', 'T', 'C', 'G', 'N'] and  altbase in ['A', 'T', 'C', 'G', 'N']:
                        #    refcount=good_count[refbase]
                        #    altcount=good_count[altbase]
                        #elif altbase=='-' or (len(refbase) > 1 and len(altbase)==1):
                        #    refcount=good_readN
                        #    altcount=del_readN
                        #elif refbase=='-' or (len(altbase) > 1 and len(refbase)==1):
                        #    refcount=ins_readN
                        #    altcount=good_readN
                        if len(refbase) > 1 and len(altbase)==1:
                            refcount=del_count['N']
                            altcount=del_count['Y']
                        elif len(altbase)>1 and len(refbase)==1:
                            refcount=ins_count['N']
                            altcount=ins_count['Y']
                        elif len(altbase)>1 and len(refbase)>1:
                            refcount=mm_count['N']
                            altcount=mm_count['Y']
                        else:
                            refcount=good_count[refbase]
                            altcount=good_count[altbase]
                        print(refcount, altcount)
                        out=map(str, [id, chr, pos1based, column.n, good_count_no_overlapN, refbase, altbase, Tag1, Tag2, refcount, altcount,\
                                       good_readN, good_count['A'], good_count['T'], good_count['C'], good_count['G'],\
                                       ins_count['Y'], del_count['Y'],\
                                       lowmapq_readN, lowmapq_count['A'], lowmapq_count['T'], lowmapq_count['C'], lowmapq_count['G'],\
                                       nomapq_readN, nomapq_count['A'], nomapq_count['T'], nomapq_count['C'], nomapq_count['G'],\
                                       lowreadq_readN, lowreadq_count['A'], lowreadq_count['T'], lowreadq_count['C'], lowreadq_count['G'],\
                                       lowbaseq_readN, lowbaseq_count['A'], lowbaseq_count['T'], lowbaseq_count['C'], lowbaseq_count['G'] ])
                        fh.write('\t'.join(out) + '\n') 
                 #return good_count_no_overlap, good_fcount, good_rcount, ins_count, ins_fcount, ins_rcount, del_count, del_fcount, del_rcount, mm_count, mm_fcount, mm_rcount
    if write_output:fh.close()
    if print_reads: fh2.close()
    #return {'A':0, 'T':0, 'C':0, 'G':0, 'N':0},{'A':0, 'T':0, 'C':0, 'G':0, 'N':0},{'A':0, 'T':0, 'C':0, 'G':0, 'N':0}, \
    #       {'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0},{'Y':0, 'N':0}
    if print_reads and barcode_reads is not None:
        add_barcode(output2, barcode_reads)
    elif barcode_reads is not None:
        print('[ERROR] Need to add print_reads to make barcode_reads work')

def add_barcode(read_fn, barcode_reads, index_len=8):
    read_dict=dict()
    with open(read_fn, 'r') as fin:
        for line in fin:
            qname, pos1based, readbase=line.strip().split('\t')
            if qname not in read_dict: read_dict[qname]=dict()
            if pos1based not in read_dict[qname]: 
                 read_dict[qname][pos1based]=dict()
                 read_dict[qname][pos1based][1]=readbase
            else:
                 read_n=max(list(read_dict[qname][pos1based].keys()))+1
                 read_dict[qname][pos1based][read_n]=readbase
                 #if readbase != read_dict[qname][pos1based]:
                 #    print('[ERROR] conflict base in read {0}: {1} vs. {2}'.format(qname, read_dict[qname][pos1based], readbase))
    #print(read_dict)
    with open(read_fn+'.with_barcode.txt', 'w') as fout, pysam.FastxFile(barcode_reads) as readin:
         headers=['read','pos','base','paired_read_n','full_index','index','UMI']
         fout.write('\t'.join(list(map(str, headers)))+'\n')
         for entry in readin:
            if len(read_dict.keys())==0:break
            if entry.name in read_dict:
                for pos in read_dict[entry.name]:
                    for read_n in read_dict[entry.name][pos]:
                        index=entry.sequence[:index_len]
                        UMI=entry.sequence[index_len:]
                        out=[entry.name, pos, read_dict[entry.name][pos][read_n],read_n, entry.sequence, index, UMI]
                        #print(out) 
                        fout.write('\t'.join(list(map(str, out)))+'\n')
                read_dict.pop(entry.name,None)
class REFERENCE:
  def __init__(self, genome_build):
    self.gb=genome_build
    self.AUTO_DIR="/rgs01/project_space/cab/automapper/common/yhui/Ti-Cheng"
    self.REF_DIR_H=self.AUTO_DIR + '/REF/Homo_sapiens/NCBI'
    self.REF_DIR_M=self.AUTO_DIR + "/REF/Mus_musculus/Gencode"
  def get(self):
    return self.genome_reference(), self.mappability(), self.gc(), self.repeatmasker(), self.CpG(), self.TFBS(), self.twobit()
  def genome_reference(self):
    REFERENCE={'hg19': self.REF_DIR_H+"/hg19_test/GRCh37-lite.fa",
               'hg19_withchr': self.REF_DIR_H+"/hg19_test/GRCh37-lite_wchr.fa",
               'hg38': self.REF_DIR_H+"/hg38_test/GRCh38_no_alt.fa",
               'mm9' : self.REF_DIR_M+"/M1/bwa-index/0.7.17-r1188/NCBIM37.genome.fa",
               'mm9_compbio': self.REF_DIR_M+"/M1/bwa-index/0.7.17-r1188/MGSCv37.fa",
               'mm10': self.REF_DIR_M+"/M22/bwa-index/0.7.17-r1188/GRCm38.primary_assembly.genome.fa"}
    if self.gb in REFERENCE:
        logging.info("[REF] "+REFERENCE[self.gb])
        return REFERENCE[self.gb]
    else:
        sys.exit('Unknown reference ... exit')

  def twobit(self):
    TWOBIT={'hg19':self.REF_DIR_H+"/hg19_test/GRCh37-lite.2bit",
            'hg19_withchr':self.REF_DIR_H+"/hg19_test/GRCh37-lite_wchr.2bit",
            'hg38':self.REF_DIR_H+"/hg38_test/GRCh38_no_alt.2bit",
            'mm9': self.REF_DIR_M+"/M1/bwa-index/0.7.17-r1188/NCBIM37.genome.2bit",
            'mm9_compbio':self.REF_DIR_M+"/M1/bwa-index/0.7.17-r1188/MGSCv37.2bit",
            'mm10': self.REF_DIR_M+"/M22/bwa-index/0.7.17-r1188/GRCm38.2bit"}
    if self.gb in TWOBIT:
        logging.info("[2BIT] "+TWOBIT[self.gb])
        return TWOBIT[self.gb]
    else:
        logging.info("[No 2BIT] no internal 2bit identified for" + self.gb )
        return None

  def mappability(self):
    MAP={'hg38': self.REF_DIR_H+"/hg38_test/hg38_100.mappability.bw",
         'hg19': self.REF_DIR_H+"/hg19_test/wgEncodeCrgMapabilityAlign100mer.bw",
         'hg19_withchr': self.REF_DIR_H+"/hg19_test/wgEncodeCrgMapabilityAlign100mer.bw"
        }
    if self.gb in MAP:
        logging.debug("[MAP] "+MAP[self.gb])
        return MAP[self.gb]
    else:
        logging.error("[no MAP] no mappability track for annotation for "+self.gb )
        return None

  def gc(self):
    GC={'hg38': self.REF_DIR_H+"/hg38_test/hg38.gc5Base.bw",
        'hg19': self.REF_DIR_H+"/hg19_test/hg19.gc5Base.bw",
        'hg19_withchr': self.REF_DIR_H+"/hg19_test/hg19.gc5Base.bw"
        }
    if self.gb in GC:
        logging.debug("[GC] "+GC[self.gb])
        return GC[self.gb]
    else:
        logging.error("[no GC] no GC track for annotation for "+self.gb )
        return None

  def repeatmasker(self):
    REPEAT={'hg38': self.REF_DIR_H+"/hg38_test/hg38_rmsk.sorted.bed.gz",
            'hg19': self.REF_DIR_H+"/hg19_test/hg19_rmsk.sorted.bed.gz",
            'hg19_withchr': self.REF_DIR_H+"/hg19_test/hg19_rmsk.sorted.bed.gz"
           }
    if self.gb in REPEAT:
        logging.debug("[REPEAT] "+REPEAT[self.gb])
        return REPEAT[self.gb]
    else:
        logging.error("[no REPEAT] no REPEATMASKER track for annotation for "+self.gb )
        return None

  def CpG(self):
    CpG={'hg38': self.REF_DIR_H+"/hg38_test/hg38_CpG_island.sorted.bed.gz",
         'hg19': self.REF_DIR_H+"/hg19_test/hg19_CpG_island.sorted.bed.gz",
         'hg19_withchr': self.REF_DIR_H+"/hg19_test/hg19_CpG_island.sorted.bed.gz"
        }
    if self.gb in CpG:
        logging.debug("[CpG] "+CpG[self.gb])
        return CpG[self.gb]
    else:
        logging.error("[no CpG] no CpG island track for annotation for "+self.gb )
        return None
  def TFBS(self):
    TFBS={'hg38': self.REF_DIR_H+"/hg38_test/hg38_TFBS_conserved.sorted.bed.gz",
          'hg19': self.REF_DIR_H+"/hg19_test/hg19_TFBS_conserved.sorted.bed.gz",
          'hg19_withchr': self.REF_DIR_H+"/hg19_test/hg19_TFBS_conserved.sorted.bed.gz"
        }
    if self.gb in TFBS:
        logging.debug("[TFBS] "+TFBS[self.gb])
        return TFBS[self.gb]
    else:
        logging.error("[no TFBS] no TFBS track for annotation for "+self.gb )
        return None
  def exon(self):
    EXON={'hg38':self.REF_DIR_H+"/hg38_test/hg38.exon.sorted.bed.gz",
          'hg19':self.REF_DIR_H+"/hg19_test/hg19.exon.sorted.bed.gz",
          'hg19_withchr':self.REF_DIR_H+"/hg19_test/hg19.exon.sorted.bed.gz"
         }
    if self.gb in EXON:
        logging.debug("[EXON] "+EXON[self.gb])
        return EXON[self.gb]
    else:
        logging.error("[no EXON] no EXON track for annotation for "+self.gb )
        return None
  def gene(self):
    GENE={'hg38':self.REF_DIR_H+"/hg38_test/hg38.gene.sorted.bed.gz",
          'hg19':self.REF_DIR_H+"/hg19_test/hg19.gene.sorted.bed.gz",
          'hg19_withchr':self.REF_DIR_H+"/hg19_test/hg19.gene.sorted.bed.gz"}
    if self.gb in GENE:
        logging.debug("[GENE] "+GENE[self.gb])
        return GENE[self.gb]
    else:
        logging.error("[no GENE] no GENE track for annotation for "+self.gb )
        return None
  def cytoband(self):
    CYTOBAND={'hg38':self.REF_DIR_H+"/hg38_test/hg38.cytoband.bed.gz",
              'hg19':self.REF_DIR_H+"/hg19_test/hg19.cytoband.bed.gz",
              'hg19_withchr':self.REF_DIR_H+"/hg19_test/hg19.cytoband.bed.gz"}
    if self.gb in CYTOBAND:
        logging.debug("[CYTOBAND] "+CYTOBAND[self.gb])
        return CYTOBAND[self.gb]
    else:
        logging.error("[no CYTOBAND] no CYTOBAND track for annotation for "+self.gb )
        return None



def readNSum(df):
    '''
        calcualte sum from a dictionary
    '''
    sum=0
    #for k, v in df.iteritems(): ### python 2.7
    for k, v in df.items():
        sum+=v
    return sum

def mean(list_of_quals):
    """
    Naive function for determining 'read quality'
    :param list_of_quals: list of numeric quality values
    :return: mean of the quality values
    """
    return float(sum(list_of_quals)) / len(list_of_quals)

def main():
    parser = argparse.ArgumentParser(description='single cell analses')
    parser.add_argument('-id' , '--id', help='ID', required=True)
    parser.add_argument('-bam' , '--bam', help='BAM file', required=True)
    parser.add_argument('-var' , '--var', help='Variation file with 4 columns [chr,position,ref,alt] separated by tab', required=True)
    parser.add_argument('-out' , '--out', help='Output', required=True)
    parser.add_argument('-gb' , '--genome_build', help='Genome version[hg19/hg38/mm10/mm9]', required=True)
    #parser.add_argument('-tag' , '--tag', help='Tag to collpase reads, e.g. UB (error corrected UMI, default), CB (error corrected barcode), UR (UMI), CR(barcode)', default='UB')
    parser.add_argument('-bq' , '--minimum_base_quality', help='minimum per base quality cutoff (default: 20)', default=20, type=int)
    parser.add_argument('-rq' , '--minimum_read_quality', help='minimum mean read quality cutoff (default: 20)', default=20, type=int)
    parser.add_argument('-mq' , '--minimum_mapping_quality', help='minimum mapping quality cutoff (default: 10)', default=10, type=int) ### 255 is the default for cellranger
    parser.add_argument('-dedup' , '--de_duplicate', help='de-deuplicate reads', action='store_true')
    parser.add_argument('-preads' , '--print_reads', help='whether to print reads', action='store_true')
    parser.add_argument('-preads_b' , '--print_reads_bases', help='print only the reads harboring these bases, separated by comma')
    parser.add_argument('-skip_indel' , '--skip_indel', help='skip indel, the counting would fail for the xenocp BAMs', action='store_true')
    parser.add_argument('-bc_reads' , '--barcode_reads', help='barcode reads', required=False) 

    args=parser.parse_args() 
    if not args.genome_build:
            sys.exit('Please provide genome build or reference ... exit')
    else:
            args.reference, args.mappability, args.gc, args.repeatmasker, args.cpg, args.tfbs, twobit=REFERENCE(args.genome_build).get()

 
    var_read_depth(args.var, args.bam, args.out, args.id, args.reference, \
            minimum_base_quality = args.minimum_base_quality, minimum_mapq = args.minimum_mapping_quality, minimum_read_quality = args.minimum_read_quality, dedup=args.de_duplicate, print_reads=args.print_reads, print_base=args.print_reads_bases, stepper='all',skip_indel=args.skip_indel, barcode_reads=args.barcode_reads)


if __name__ == "__main__":
    main()                                                       
