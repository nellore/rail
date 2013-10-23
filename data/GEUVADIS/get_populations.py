# get geuvadis population data from ENA + HapMap websites
# AF 23 Oct 2013

#### download ENA accesson numbers 
import urllib
import urllib2

ena_info = urllib2.urlopen("http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22secondary_study_accession=%22ERP001942%22%22&result=read_run&limit=667&length=667&offset=0&display=report&fields=study_accession,secondary_study_accession,sample_accession_list,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,col_taxonomy")
dat = ena_info.readlines()
table_rows = [x.split('\t') for x in dat]
sample_ftps = [x[10] for x in table_rows]
len(set(sample_ftps)) #668, too many
numinfront = len("ftp.sra.ebi.ac.uk/vol1/ERA169/ERA169774/fastq/")
hapmapids = [x[numinfront:][:7] for x in sample_ftps]
len(set(hapmapids)) #464 unique hapmap samples (-1 for col name)

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
with open(out_table, 'w') as f:
    f.write('run_id\thapmap_id\tpopulation\n')
    tr = table_rows[1:]
    pos = 0
    for h in hapmapids[1:]:
        f.write(tr[pos][4]+'\t'+h+'\t'+lookup[h]+'\n')
        pos += 1





