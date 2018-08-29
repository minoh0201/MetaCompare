#!/usr/bin/python3

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-' and argv[0][1] != 'h':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def intersection(list1, list2):
    return list(set(list1) & set(list2))

if __name__ == '__main__':
    from sys import argv
    from sys import exit
    import os
    import subprocess
    import sys
    from Bio import SeqIO
    import pandas as pd
    import math

    myargs = getopts(argv)
    if '-h' in myargs or len(myargs) == 0:  # Example usage.
        print('\nUsage: ./metacmp.py -c filename1.fa -g filename2.fa [-t 64] \n')
        print('\t-c: Specify FASTA file containing assembled contigs (resulted from IDBA-UD).')
        print('\t-g: Specify FASTA file containing predicted genes (derived from prodigal).')
        print('\t-t: Specify the number of threads will be used in executing blast (default: 64).')
		print('\t-v: Printing important values in calculation On/Off [1: On (default), 0: Off]')
        print()
        exit()

    if not '-t' in myargs:
        myargs['-t'] = '64'
	
	if not '-v' in myargs:
        myargs['-v'] = '1'

    sample_name = myargs['-c'].split('.')[0]

    acc_name = sample_name + "_ACLAME.txt"
    if not os.path.exists(os.getcwd()+"/"+acc_name):
        print('Running blastn on ACLAME')
        subprocess.call(["blastn", "-db", "./BlastDB/aclame", "-query", myargs['-c'], \
                         "-out", acc_name, "-outfmt", "6", \
                         "-num_threads", myargs['-t'], "-evalue", "1e-10"])
    else:
        print('Skipping: Blast output is already exist for ACLAME')

    card_name = sample_name + "_CARD.txt"
    if not os.path.exists(os.getcwd()+"/"+card_name):
        print('Running blastx on CARD')
        subprocess.call(["blastx", "-db", "./BlastDB/CARD_PROT", "-query", myargs['-g'], \
                         "-out", card_name, "-outfmt", "6", \
                         "-num_threads", myargs['-t'], "-evalue", "1e-10"])
    else:
        print('Skipping: Blast output is already exist for CARD')

    patric_name = sample_name + "_PATRIC.txt"
    if not os.path.exists(os.getcwd()+"/"+patric_name):
        print('Running blastn on PATRIC')
        subprocess.call(["blastn", "-db", "./BlastDB/PATRIC", "-query", myargs['-c'], \
                         "-out", patric_name, "-outfmt", "6", \
                         "-num_threads", myargs['-t'], "-evalue", "1e-10"])
    else:
        print('Skipping: Blast output is already exist for PATRIC')

    print('Reading files..')

    #Open Fasta file
    records = list(SeqIO.parse(myargs['-c'], "fasta"))
    nContigs = len(records)
    #print(records[0].id)
    #print(records[0])

    # blast output column name
    # ## Fields: 
    #0: query acc.ver
    #1: subject acc.ver
    #2: % identity
    #3: alignment length
    #4: mismatches
    #5: gap opens
    #6: q. start
    #7: q. end
    #8: s. start
    #9: s. end
    #10: evalue
    #11: bit score


    print('Computing resistome risk score..')
 
    # Open blast output for ACLAME
    if not os.path.getsize(os.getcwd()+"/"+acc_name) > 0:
        #file is empty
        print('Warning: '+ acc_name+ ' is empty.')
        MGE_contigs = []
    else:
        ac = pd.read_csv(acc_name, sep='\t', header=None)
        ac.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
        # filter out contigs identity under 60
        ac_iden60 = ac[ac.identity > 60]

        # filter out contigs having less than 90% coverage of the reference
        if not os.path.exists(os.getcwd()+"/Len_aclame.txt"):
            print('Len_aclame.txt file does not exists.')
            sys.exit()
        else:
            acleng = pd.read_csv("Len_aclame.txt", sep='\t', header=None)
            acleng.columns = ['sub_id', 'ref_gene_leng']
            ac_merged = pd.merge(ac_iden60, acleng, how = 'left', on = 'sub_id')
            ac_filtered = ac_merged[ac_merged.alignLen > (ac_merged.ref_gene_leng * 0.9)]

            # Note: 'ac_filtered' data frame contains name and position of MGEs in contigs
            MGE_contigs = ac_filtered.id.unique()

    # Open blast output for PATRIC
    if not os.path.getsize(os.getcwd()+"/"+patric_name) > 0:
        #file is empty
        print('Warning: '+ patric_name+ ' is empty.')
        PAT_contigs = []
    else:
        pa = pd.read_csv(patric_name, sep='\t', header=None)
        pa.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
        # filter out contigs identity under 60
        pa_iden60 = pa[pa.identity > 60]

        # filter out contigs alignment length under 150
        pa_alen = pa_iden60[pa_iden60.alignLen > 150]

        # Note: 'pa_alen' data frame contains name and position of Pathogens in contigs
        PAT_contigs = pa_alen.id.unique()

    # Open blast output for CARD
    if not os.path.getsize(os.getcwd()+"/"+ card_name) > 0:
        #file is empty
        print('Warning: '+ card_name+ ' is empty.')
        CARD_contigs = []
    else:
        ca = pd.read_csv(card_name, sep='\t', header=None)
        ca.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
        # filter out contigs identity under 60
        ca_iden60 = ca[ca.identity > 60]

        # filter out contigs alignment length under 25
        ca_alen = ca_iden60[ca_iden60.alignLen > 25]

        # Note: 'ca_alen' data frame contains name and position of Antibiotic Resistence Genes in contigs
        CARD_genes = ca_alen.id.unique()
        ARG_contigs = []
        for gene in ca_alen.id:
            contig_name = gene.split("_")
            ARG_contigs.append(contig_name[0]+"_"+contig_name[1])
        ARG = set(ARG_contigs)

        # computing risk score
        ARG_MGE = intersection(ARG_contigs, MGE_contigs)
        ARG_MGE_PAT = intersection(ARG_MGE, PAT_contigs)
        nARG = len(ARG)
        nARG_MGE = len(ARG_MGE)
        nARG_MGE_PAT = len(ARG_MGE_PAT)

        fARG = float(nARG)/nContigs
        fARG_MGE = float(nARG_MGE)/nContigs
        fARG_MGE_PAT = float(nARG_MGE_PAT)/nContigs

        distance = math.sqrt((0.01 - fARG)**2 + (0.01 - fARG_MGE)**2 + (0.01 - fARG_MGE_PAT)**2)

        score = 1.0 / ( (2 + math.log10(distance))**2 )

		# other stat

		nMGE = len(MGE_contigs)
		nPAT = len(PAT_contigs)
		fMGE = float(nMGE)/nContigs
		fPAT = float(nPAT)/nContigs

        print("Resistome risk score: " + str(score))

        if myargs['-v'] == '1':
            print("nContigs, nARG, nMGE, nPAT, nARG&MGE, nARG&MGE&PAT, nARG/nContigs, nMGE/nContigs, nPAT/nContigs, nARG&MGE/nContigs, nARG&MGE&PAT/nContigs, Risk Score\n")
            print(nContigs, nARG, nMGE, nPAT, nARG_MGE, nARG_MGE_PAT, fARG, fMGE, fPAT, fARG_MGE, fARG_MGE_PAT, score)





















