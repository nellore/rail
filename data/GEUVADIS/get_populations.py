# get geuvadis population data from ENA + HapMap websites
# AF 23 Oct 2013

#### download ENA accesson numbers 
import urllib2
ena_info = urllib2.urlopen("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERP001942&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,col_tax_id,col_scientific_name")
#ena_info = urllib2.urlopen("http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22secondary_study_accession=%22ERP001942%22%22&result=read_run&limit=667&length=667&offset=0&display=report&fields=study_accession,secondary_study_accession,sample_accession_list,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,col_taxonomy") ### OUTDATED as of 2/8/14
dat = ena_info.readlines()
table_rows = [x.split('\t') for x in dat][1:] #strip off header row
col_labels = dat[0].split('\t')
ftpind = col_labels.index('submitted_ftp') #11
sample_ftps = [x[ftpind] for x in table_rows]
len(set(sample_ftps)) #667
numinfront = len("ftp.sra.ebi.ac.uk/vol1/ERA169/ERA169774/fastq/")
hapmapids = [x[numinfront:][:7] for x in sample_ftps]
len(set(hapmapids)) #464 unique hapmap samples 

## note: from paper supplement:
## 5 samples sequenced in replicate at each of 7 labs (so 8x total)
## 168 samples sequenced in replicate at 2/3 coverage
## sooo, we have 465 base samples + 35 + 168 replicates = 668
## one sample is missing from ENA.

sample_info = urllib2.urlopen("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/sequence.index")
popdat = sample_info.readlines()
pdrows = [x.split('\t') for x in popdat]
idcol = pdrows[0].index('SAMPLE_NAME')
popcol = pdrows[0].index('POPULATION')
lookup = dict()
for entry in pdrows[1:]:
    lookup[entry[idcol]] = entry[popcol]

out_table = 'pop_data.txt'
errind = col_labels.index('run_accession')
with open(out_table, 'w') as f:
    f.write('run_id\thapmap_id\tpopulation\n')
    pos = 0
    for h in hapmapids:
        f.write(table_rows[pos][errind]+'\t'+h+'\t'+lookup[h]+'\n')
        pos += 1

# running a diff on the output of this and the old copy of pop_data showed that everything from before is still right - AWESOME

# now add identifiers to this table (to look up in the table from the geuvadis authors)
ids = [x.split(';')[0].split('/')[-1][:-11] for x in sample_ftps]
with open(out_table, 'w') as f:
    f.write('sample_id\trun_id\thapmap_id\tpopulation\n')
    pos = 0
    for h in hapmapids:
        f.write(ids[pos]+'\t'+table_rows[pos][errind]+'\t'+h+'\t'+lookup[h]+'\n')
        pos += 1








