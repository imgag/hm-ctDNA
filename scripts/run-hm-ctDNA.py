import logging
import os
import warnings
import pandas as pd
import subprocess as sp
import sys
import shlex
import datetime
import yaml

# basic configuration
root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

# set basic pandas options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)
env = os.environ.copy()
env["PATH"] = "/path/to/methyldackel/bin:" + env["PATH"]

##
steps = ['bam','md','bedpe','combine','sum1','stat']

##config file
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
config_path = os.path.join(base_dir, "config", "config.yaml")

with open(config_path) as f:
    config = yaml.safe_load(f)

opt_folder = config['ngs_bits_folder']
project_folder = config["project_folder"]
bam_folder = config["bams"]
tmp_bed = config['bed']
ref_fa = config['ref_fa']
out_folder = config["out_folder"]

# create out_folder
os.makedirs(out_folder, exist_ok=True)

# read sample file
sample_file = config['meta_data']
df_samples = pd.read_csv(sample_file, sep="\t", index_col=False)

cpgs = []
df_bed = pd.read_csv(tmp_bed, sep="\t", header=None, index_col=False)

for i, r in df_bed.iterrows():
    entry = r[3]  # e.g. cluster002:Pan-BRCA_Hyper or cpg123456:Pan-BRCA_Hyper

    try:
        cluster_id, cluster_group = entry.split(':')
    except ValueError:
        logging.warning(f"Unexpected format: {entry}")
        continue

    cpgs.append(cluster_id)
# optional: remove duplicates
cpgs = list(set(cpgs))

logging.info(f"Loaded {len(cpgs)} clusters")

## prep 2 - check if all needed files are available, copy, filter, sort and index bam file
rows = []
count = 1
logging.info('Check if all input files are available and copy bam file...')
print(df_samples.columns.tolist())
for i,r in df_samples.iterrows():
    sample = r['#Sample-ID']
    logging.info(' Processing sample {} ({}/{})'.format(sample,count,len(df_samples.index)))
    count += 1
    
    umi = False
    if r['UMI'] == 'yes':
        umi = True

    bam = None
    if umi:
        bam = '{}/{}.sorted.bam'.format(bam_folder, sample)
        if not os.path.isfile(bam):
            logging.warning(' Could not find file {}. Skipping.'.format(bam))
            continue
    else:
        bam = '{}/{}.sorted.bam'.format(bam_folder, sample)

    tmp_bam0 = '{}/{}-tmp0.bam'.format(out_folder, sample)
    tmp_bam1 = '{}/{}-tmp1.bam'.format(out_folder, sample)
    tmp_bam2 = '{}/{}-tmp2.bam'.format(out_folder, sample)
    tmp_bam3 = '{}/{}-tmp3.bam'.format(out_folder, sample)
    
    # 2. filter bam file for target regions and extract paired reads
    if 'bam' in steps:
        # check if index file is available or generate one
        if not os.path.isfile('{}.bai'.format(bam)):
            logging.info(' Indexing file {}.'.format(bam))
            cmd = 'samtools index {}'.format(bam)
            # print(cmd)
            sp.check_output(shlex.split(cmd), stderr=sp.STDOUT)

        # extract pairs overlapping target regions
        # h - with header, b - bam format, P - fetch pairs, L - target bed file
        cmd = 'samtools view -hbP -L {} -o {} {}'.format(tmp_bed,tmp_bam0,bam)
        sp.run(shlex.split(cmd))
        
        ###filter_fraglength
        tmp_bam0_filt = '{}/{}-tmp0_filt.bam'.format(out_folder, sample)
        cmd = "samtools view -h {} | awk '($9 >= 125 && $9 <= 180) || ($9 <= -125 && $9 >= -180) || $1 ~ /^@/' | samtools view -b -o {}".format(tmp_bam0, tmp_bam0_filt)
        #print(cmd)
        sp.run(cmd, shell=True)

        ## index bam file
        cmd = 'samtools index {}'.format(tmp_bam0)
        sp.check_output(shlex.split(cmd),stderr=sp.STDOUT)

        # sort bam by name in preparation for overlap clipping
        cmd = 'samtools sort -n -o {} {}'.format(tmp_bam1,tmp_bam0_filt)
        #print(cmd)
        sp.check_output(shlex.split(cmd),stderr=sp.STDOUT)

        # overlap clipping
        cmd = '{}/ngs-bits-hg38-2023_02-1-g2a03d5a8/BamClipOverlap -in {} -out {}'.format(opt_folder,tmp_bam1,tmp_bam2)
        #print(cmd)
        sp.run(shlex.split(cmd))

        # sort bam
        cmd = 'samtools sort -o {} {}'.format(tmp_bam3,tmp_bam2)
        #print(cmd)
        sp.check_output(shlex.split(cmd),stderr=sp.STDOUT)

        # index bam file
        cmd = 'samtools index {}'.format(tmp_bam3)
        #print(cmd)
        sp.check_output(shlex.split(cmd),stderr=sp.STDOUT)

    else:
        if not os.path.isfile(tmp_bam1):
            logging.warning(' File {} missing. Skipping.'.format(tmp_bam1))
            continue

    rows.append(r.to_list() + [tmp_bam3,tmp_bam2])

columns = df_samples.columns.to_list()
columns.append('BAM')
columns.append('BAM-NAME-SORTED')
df_samples = pd.DataFrame(data=rows, columns=columns)

## analysis 1 run MethylDackel on prepared bam
if 'md' in steps:
    logging.info('Per read analysis with MethylDackel...')

    count = 1
    for i,r in df_samples.iterrows():
        sample = r['#Sample-ID']
        # 1. run MethylDackel perRead
        logging.info(' Analysing sample {} ({}/{})...'.format(sample,count,len(df_samples.index)))
        count += 1
        bam = r['BAM']
        out_file = f'{out_folder}/{sample}-tmp2.tsv'
        cmd = f"{config['Methyldackel_path']} perRead {ref_fa} {bam} -o {out_file}"

        #print(cmd)
        sp.check_output(shlex.split(cmd),stderr=sp.STDOUT, env=env)

## analysis 2 generate bedpe files
if 'bedpe' in steps:
    logging.info('Bam to bedpe...')

    count = 1
    for i,r in df_samples.iterrows():
        sample = r['#Sample-ID']
        bam = r['BAM-NAME-SORTED']
        logging.info(' Processing sample {} ({}/{})'.format(sample,count,len(df_samples.index)))
        count += 1

        # 3. convert bam to bedpe and filter for target bed file
        cmd = []
        cmd.append('{}/bedtools2/bin/pairToBed -bedpe -abam {} -b {}'.format(opt_folder,bam,tmp_bed))
        cmd.append('sort -u')
        tmp_out = open('{}/{}-tmp4.bedpe'.format(out_folder,sample), 'w')
        #print('{} > {}'.format(' | '.join(cmd),tmp_out.name))
        ps1 = sp.Popen(shlex.split(cmd[0]),stdout=sp.PIPE)
        ps2 = sp.Popen(shlex.split(cmd[1]),stdin=ps1.stdout,stdout=tmp_out)
        #ps1.wait()
        ps2.wait()
        tmp_out.close()

## analysis 3 combine bedpe and MethylDackel results
if 'combine' in steps:
    logging.info('Combine bedpe and MethylDackel...')
    
    count = 1
    for i, r in df_samples.iterrows():
        sample = r['#Sample-ID']
        logging.info(' Processing sample {} ({}/{})'.format(sample,count,len(df_samples.index)))
        count += 1
        # 4. combine MethylDackel output and bedpe file
        
        # sort MethylDackel output with the read name column; in this case column 1
        cmd = 'sort -k1 {}/{}-tmp2.tsv -o {}/{}-tmp2_sorted.tsv'.format(out_folder,sample,out_folder,sample)
        #print(cmd)
        sp.run(shlex.split(cmd),stderr=sp.STDOUT)
        # sort pairToBed output with sample with the read name; in this case column 7
        cmd = 'sort -k7 {}/{}-tmp4.bedpe -o {}/{}-tmp4_sorted.bedpe'.format(out_folder,sample,out_folder,sample)
        #print(cmd)
        sp.run(shlex.split(cmd),stderr=sp.STDOUT)
        # resort pairToBed file to match the TsvMerge needs (both files with read name columns as first columns)
        cmd = '{}/ngs-bits-hg38-2023_02-1-g2a03d5a8/TsvSlice -numeric -cols 7,1,2,3,4,5,6,7,8,9,14 -in {}/{}-tmp4_sorted.bedpe -out {}/{}-tmp4_sliced.bedpe'.format(opt_folder,out_folder,sample,out_folder,sample,out_folder,sample)
        #print(cmd)
        sp.run(shlex.split(cmd))

        # 5. remove duplicate regions
        file1 = open('{}/{}-tmp4_sliced.bedpe'.format(out_folder,sample), 'r')
        file2 = open('{}/{}-tmp4_dedup.bedpe'.format(out_folder,sample), 'w')
        prev_cols = None
        cur_cols = None
        while True:
            # get next line from file
            cur_line = file1.readline(-1)
            # if line is empty end of file is reached
            if not cur_line:
                break

            cur_cols = cur_line.replace('\n', '').split('\t')
            if not prev_cols:
                prev_cols = cur_cols
                continue
            if (cur_cols[0:9] == prev_cols[0:9]):
                prev_cols[10] = '{};{}'.format(prev_cols[10],cur_cols[10])
                continue

            file2.write("{}\n".format('\t'.join(prev_cols)))
            prev_cols = cur_cols

        file2.write("{}\n".format('\t'.join(prev_cols)))
        file1.close()
        file2.close()

        # merge pairToBed file and MethylDackel file using read names
        cmd = '{}/ngs-bits-hg38-2023_02-1-g2a03d5a8/TsvMerge -numeric -cols 1 -in {}/{}-tmp4_dedup.bedpe {}/{}-tmp2_sorted.tsv -out {}/{}-tmp5.tsv'.format(opt_folder,out_folder,sample,out_folder,sample,out_folder,sample)
        #print(cmd)
        sp.run(shlex.split(cmd))

# analysis 3 combine MethylDackel and BedPe output

error1 = open('{}/{}-tmp6_errors.txt'.format(out_folder, sample), 'w')
if 'sum1' in steps:
    logging.info('Summary of bedpe files...')
    count = 1
    for i,r in df_samples.iterrows():
        sample = r['#Sample-ID']
        logging.info(' Processing sample {}'.format(sample,count,len(df_samples.index)))
        count += 1

        counts = {}
        #for ob in oncogenic_brca:
        #    counts[ob] = 0
        for cpg in cpgs:
            counts[cpg] = 0

        # 7. resort cols of bedpe file and filter for methylation regions
        bedpe_lines = []
        file1 = open('{}/{}-tmp5.tsv'.format(out_folder,sample), 'r')
        file2 = open('{}/{}-tmp6.bedpe'.format(out_folder,sample), 'w')

        count = 0
        count_oncogenic = 0
        # write header for file 2
        header = ['chr_1','start_1','end_1','chr_2','start_2','end_2','','strand','read_id','cpg_id','percm_cpg','count_cpg']
        file2.write("#{}\n".format('\t'.join(header)))
        while True:
            # Get next line from file
            line = file1.readline(-1)

            # if line is empty end of file is reached
            if not line:
                break
            # skip comments (ngs-bits adds a header
            if line.startswith('#') or line.startswith('##'):
                continue

            columns = line.replace('\n', '').split('\t')
            #print(columns)

            # empty first part, this means result of MethylDackel is there but not Bedtools
            if not columns[1]:
                error1.write("Empty first part\n{}\n".format('\t'.join([str(x) for x in columns])))
                continue
            # empty second part, this means result of Bedtoosl is there but not MethylDackel
            if not columns[11]:
                error1.write("Empty second part\n{}\n".format('\t'.join([str(x) for x in columns])))
                continue

            myorder = []
            (idx_cpgs_1, idx_percm_1, idx_cpgs_2, idx_percm_2) = [None,None,None,None]
            if len(columns)==19:
                (idx_cpgs_1,idx_percm_1,idx_cpgs_2,idx_percm_2) = [14,13,18,17]
            else:
                logging.error('Expected number of columns, found {}. {}'.format(len(columns),' '.join(columns)))
                error1.write("Format {} columns: {}\n".format(len(columns),'\t'.join([str(x) for x in columns])))
                continue

            # calculate methylation scores and filter columns
            count_cpgs = int(columns[idx_cpgs_1]) + int(columns[idx_cpgs_2])
            percm_cpgs = 0
            if count_cpgs > 0:
                countm_cpgs = round(float(columns[idx_percm_1])/100*float(columns[idx_cpgs_1]),0) + round(float(columns[idx_percm_2])/100*float(columns[idx_cpgs_2]),0)
                percm_cpgs = countm_cpgs/count_cpgs

            myorder = [1,2,3,4,5,6,8,9,0,10]
            columns = [columns[i] for i in myorder]
            columns += [percm_cpgs,count_cpgs]
            columns = [str(x) for x in columns]
            file2.write("{}\n".format('\t'.join(columns)))

        file1.close()
        file2.close()

if 'stat' in steps:
    
    parameters = [[6, 0.00, 0.75]]
    #parameters = [[4, 0.00, 0.75], [6, 0.00, 0.75], [8, 0.00, 0.75]]

    for p in parameters:
        logging.info('Statistics ({}-{}-{})...'.format(p[0], p[1], p[2]))

        min_cpgs = p[0]
        min_hyper_meth_frac = p[2]

    # initialize output dataframe
        df = pd.DataFrame()
        df['#CPGs'] = cpgs

        for sample in df_samples['#Sample-ID'].to_list():
            logging.info(f'sample {sample}')

        # initialize counts
            counts_true = {c: 0 for c in cpgs}
            counts_false = {c: 0 for c in cpgs}

            with open(f'{out_folder}/{sample}-tmp6.bedpe', 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    columns = line.strip().split('\t')

                    percm_cpg = float(columns[10])
                    count_cpg = int(columns[11])

                # filter low CpG fragments
                    if count_cpg < min_cpgs:
                        continue

                # split multiple CpGs
                    tmp_cpgs = columns[9].split(';')

                    for entry in tmp_cpgs:
                        tmp_cpg = entry.split(':')[0].strip()

                        if tmp_cpg not in counts_true:
                            continue

                    # classify
                        if percm_cpg >= min_hyper_meth_frac:
                            counts_true[tmp_cpg] += 1
                        else:
                            counts_false[tmp_cpg] += 1

                # add to dataframe
            df[f'{sample}_oncogenic'] = df['#CPGs'].map(counts_true).fillna(0).astype(int)
            df[f'{sample}_nononcogenic'] = df['#CPGs'].map(counts_false).fillna(0).astype(int)

            # save
        out_file = f'{out_folder}/summary_all-samples-{min_cpgs}-{int(p[1]*100)}-{int(p[2]*100)}.tsv'
        df.to_csv(out_file, sep='\t', index=False)

    logging.info(f'Wrote {out_file}')

    

error1.close()
logging.info('Done.')
