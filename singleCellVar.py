#!/usr/bin/env python

import re, argparse

def summary(args, CB_ls=None, header=True, excludeNoCB=True, excludeNoUB=True, sep="\t", group_by=['ID', 'chr', 'position', 'ref', 'alt', 'CB']):
    infile=open(args.input, "r")
    rown=0
    colidx={}
    allcols=['ID','chr','position','total_n','QC_n','ref','alt','CB','UB','refcount','altcount', 'good_readN', 'ins_readN', 'del_readN', 'lowmapq_readN', 'nomapq_readN', 'lowreadq_readN', 'lowbaseq_readN']
    allcolidx={}
    summary={}

    for line in infile:
        fs=line.strip().split(sep)
        rown+=1
        grouped=[]
        allcolumns=[]
        if header and rown==1:
           hidx=0
           for i in fs:
               for g in allcols:
                   if i==g:
                       allcolidx[g]=hidx
                       allcolumns.append(g)
                       if g in group_by:
                           colidx[g]=hidx
                           grouped.append(g)
               hidx+=1
           #print '\t'.join(grouped)
           #print '\t'.join(allcolumns)
        #elif rown <= 10:
        else:
           grouped=[]
           for g in group_by:
               grouped.append(fs[colidx[g]])
           for g in allcols:
               allcolumns.append(fs[allcolidx[g]]) 

           if excludeNoCB:
               if fs[allcolidx['CB']]=='noCB': continue
           if excludeNoUB:
               if fs[allcolidx['UB']]=='noUB': continue
           if CB_ls is not None:
                compared_CB=re.sub("-.*","", fs[allcolidx['CB']])
                if compared_CB not in CB_ls:
                    ### skip the barcode not in the list
                    continue
            

           #print '\t'.join(grouped)
           #print '\t'.join(allcolumns)
           grouped_str=':'.join(grouped)
           ub=fs[allcolidx['UB']]
           if grouped_str not in summary:
                summary[grouped_str]={}
                summary[grouped_str]['nUMI']=[]
                summary[grouped_str]['ref_nUMI']=[]
                summary[grouped_str]['alt_nUMI']=[]
                if int(fs[allcolidx['refcount']]) > 0 or int(fs[allcolidx['altcount']]) > 0:
                    summary[grouped_str]['nUMI'].append(ub)
                if int(fs[allcolidx['refcount']]) > 0: 
                    summary[grouped_str]['ref_nUMI'].append(ub)
                if int(fs[allcolidx['altcount']]) > 0: 
                    summary[grouped_str]['alt_nUMI'].append(ub)
                summary[grouped_str]['ref_nReads']=int(fs[allcolidx['refcount']]) if int(fs[allcolidx['refcount']]) > 0 else 0
                summary[grouped_str]['alt_nReads']=int(fs[allcolidx['altcount']]) if int(fs[allcolidx['altcount']]) > 0 else 0
           else:
                if ub not in summary[grouped_str]['nUMI']: 
                    if int(fs[allcolidx['refcount']]) > 0 or int(fs[allcolidx['altcount']]) > 0: 
                        summary[grouped_str]['nUMI'].append(ub)
                if int(fs[allcolidx['refcount']]) > 0: 
                    if ub not in summary[grouped_str]['ref_nUMI']: summary[grouped_str]['ref_nUMI'].append(ub)
                    summary[grouped_str]['ref_nReads']+=int(fs[allcolidx['refcount']])
                if int(fs[allcolidx['altcount']]) > 0: 
                    if ub not in summary[grouped_str]['alt_nUMI']: summary[grouped_str]['alt_nUMI'].append(ub)
                    summary[grouped_str]['alt_nReads']+=int(fs[allcolidx['altcount']])
    print_summary(summary, group_by, args.output)      
           

def cells(fn, col='BC', sep="\t"):
    fh=open(fn, 'r')
    linenum=0
    headers={}
    BCs=[]
    for line in fh:
        linenum+=1
        line=line.strip()
        vals=line.split(sep)
        if linenum == 1 :
            ix=0
            for v in vals:
               headers[v]=ix
               ix+=1
        else:
            bx=vals[headers[col]]
            bx=re.sub("-.*","", bx)
            if bx not in BCs:
                BCs.append(bx)
            else:
                print "[Warning] duplicated cell barcodes:", vals[headers[col]]
    print "Number of barcodes:", len(BCs)
    return BCs


 
def print_summary(df, header, output):
    h=header
    extrah=['nUMI', 'ref_nUMI','alt_nUMI','alt_nUMI_ratio','ref_nReads','alt_nReads','alt_nReads_ratio']
    h.extend(extrah)


    #print '\t'.join(h)
    outf=open(output, 'w')
    outf.write('\t'.join(h) + '\n')
    
    for i in df:
        l=i.split(':')
        l.append(len(df[i]['nUMI']))
        refUn=len(df[i]['ref_nUMI'])
        altUn=len(df[i]['alt_nUMI'])
        l.append(refUn)
        l.append(altUn)
        if refUn+altUn > 0:
            Uratio=round(altUn/(float(refUn)+float(altUn)),5)
            l.append(Uratio)
        else:
            l.append(0)
        l.append(df[i]['ref_nReads'])
        l.append(df[i]['alt_nReads'])
        if int(df[i]['ref_nReads'])+int(df[i]['alt_nReads']) > 0:
            l.append(round(int(df[i]['alt_nReads'])/(float(df[i]['ref_nReads'])+float(df[i]['alt_nReads'])),5))
        else:
            l.append(0)
        #print '\t'.join(map(str,l))
        outf.write('\t'.join(map(str,l)) + '\n')
    outf.close()
                               
           
        
         
        

def main():
    parser = argparse.ArgumentParser(description='single cell analses, var depth')
    parser.add_argument('-i' , '--input', help='File for summarize', required=True)
    parser.add_argument('-o' , '--output', help='Output file', required=True)
    parser.add_argument('-cell', '--cell_ls', help='Limit the count to listed cells', required=False)
    args=parser.parse_args()
    cs=cells(args.cell_ls) if args.cell_ls is not None else None
    summary(args, CB_ls=cs)


if __name__ == "__main__":
    main()


