'''
Written by Heng Pan @ Weill Cornell Medicine
Feb 2018
'''

import getopt, sys
import numpy as np

try:
   opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
except getopt.GetoptError as error:
   print str(error)
   sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
       print 'USAGE: python pdrCall_from_Bismark.py [options] <sample_name> <input_dir> <output_dir>\n'
       print 'Arguments:\n'
       print '<sample_name> Input fastq names to bismark.\n'
       print '<input_dir>   Input directory where files are stored.\n'
       print '<output_dir>  Write all output files into this directory. \n'
       print 'Options:\n'
       print '-h/--help Display the help file.\n'
       sys.exit()

if len(sys.argv)-1-2*len(opts) != 3: raise EOFError('Cannot locate input file or output directors!')
sample = sys.argv[len(sys.argv)-3]
input_dir = sys.argv[len(sys.argv)-2]
output_dir = sys.argv[len(sys.argv)-1]

output_name = output_dir + '/' + 'pdr.' + sample + '.txt'
output_file = open(output_name, 'w')
line = 'chr\tstart\tstrand\tConMethReadCount\tConUMethReadCount\tDisReadCount\tNAReadCount\n'
output_file.write(line)

Meth = 0
UMeth = 1

for strand in ('OT', 'OB'):
    ss = '+' if strand == 'OT' else '-'
    input_name = input_dir + '/CpG_' + strand + '_' + sample + '_trimmed.fq_bismark.txt'
    reads = {}
    cpgs = {}

    with open(input_name) as foo:
         for i, line in enumerate(foo):
             if i % 1000000 == 0: print >> sys.stderr, 'Proccessed lines', i
             if i == 0: continue
             arr = line.strip().split('\t')

             if not (arr[2] in cpgs): cpgs[arr[2]] = {}

             if not (arr[0] in reads): reads[arr[0]] = [0, 0]
             if ((arr[1] == '+') and (arr[-1] == 'Z')): reads[arr[0]][Meth] += 1
             elif ((arr[1] == '-') and (arr[-1] == 'z')): reads[arr[0]][UMeth] += 1

             if not (int(arr[3]) in cpgs[arr[2]]): cpgs[arr[2]][int(arr[3])] = []
             cpgs[arr[2]][int(arr[3])].append(arr[0])
    foo.close()
          
    reads_id = np.asarray(reads.keys())         
    mat = np.asarray(reads.values())
    idx = np.sum(mat, axis=1) >= 4
    reads = dict(zip(reads_id[idx], mat[idx]))
    
    for read in reads:
        read_cat = ''
        if (reads[read][Meth] == 0): read_cat = 'ConUMeth'
        elif (reads[read][UMeth] == 0): read_cat = 'ConMeth'
        else: read_cat = 'Dis'
        reads[read] = read_cat
    
    for chrom in sorted(cpgs.keys()):
        for cpg in sorted(cpgs[chrom].keys()):
            ConMeth, ConUMeth, Dis, NA = 0, 0, 0, 0
            for read_id in cpgs[chrom][cpg]:
                if read_id in reads:
                   if reads[read_id] == 'ConMeth': ConMeth += 1
                   elif reads[read_id] == 'ConUMeth': ConUMeth += 1
                   elif reads[read_id] == 'Dis': Dis += 1
                else: NA += 1
            total_reads = ConMeth + ConUMeth + Dis + NA
            if total_reads >= 10:
               record = chrom + '\t' + str(cpg) + '\t' + ss + '\t' + str(ConMeth) + '\t' + str(ConUMeth) + '\t' + str(Dis) + '\t' + str(NA) + '\n'
               output_file.write(record)    
output_file.close()
