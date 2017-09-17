nordqiime study 15_Hil_CS
sequencing center - National Research Council of Canada laboratory
name as handed off from Erica: a folder in g drive called 15 Hil CS -> 15_Hil_CS
subfodler 16S
i played around with the xlsx file in Google sheets, so it looks like there are two xlsx files ?!
bacteria

16S amplification from the V4 hyper-variable region

f primer: 515F (5′-GTG CCA GCM GCC GCG GTA A-3′)
r primer: 806R (5′-GGA CTA CHV GGG TWT CTA AT-3′)
expected amplicon length: 291
file-system manipulations by mark (removing spaces in filenames, etc):
changes to mapping file:  generally adding synthetics bar-codes (necessary?  mark's python script)
making sure there's a "Description" column in the RHS most whit samples numbers (could be something else... displayed in qiime and maybe phyloseq?)
datafile like LALL-PD-KC2015-146-Control-0d-CS-16S_S1_L001_R1_001.fastq.gz

can export from multiqc:
250 nt
counts 30k -107k 

visual inspection
peak of mean qualities ~ 37... very nice and consistent

can parse out of nordqiime log:
eg:
grep "Percent combined" ~/gitstage/uderica/nordqiime/V-HMC/Fungi/nordqiime.log | sort -r  > ~/gitstage/uderica/nordqiime/V-HMC/Fungi/flashed_pct.log
%combined by flash all > 98%


extended:
lengths all close to 292 nt (multiqc table0
29k to 106k reads

try to identify unextended reads... blast against nt or refseq?

peaks of qualities ~ 38

is there anything that "doesn’t get split/picked?"

nordqiime.log messages
BEFORE CHLOROPLAST REMOVAL
Num samples: 27
Num observations: 7702
Total count: 1491582
Table density (fraction of non-zero values): 0.152

Counts/sample summary:
 Min: 27348.0
 Max: 99949.0
 Median: 54498.000
 Mean: 55243.778
 Std. dev.: 16688.446
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
191: 27348.0
187: 30582.0
202: 36991.0
207: 37267.0
197: 39670.0
177: 40259.0
172: 43116.0
147: 43285.0
198: 44505.0
201: 47439.0
206: 49254.0
186: 50695.0
188: 53473.0
203: 54498.0
196: 56743.0
178: 57915.0
208: 60043.0
192: 60189.0
183: 60663.0
173: 61205.0
182: 63024.0
193: 63969.0
146: 65173.0
176: 77683.0
171: 78456.0
148: 88188.0
181: 99949.0
Num samples: 27
Num observations: 7671
Total count: 1487431
Table density (fraction of non-zero values): 0.152

AFTER CHLOROPLAST REMOVAL
Counts/sample summary:
 Min: 27318.0
 Max: 99859.0
 Median: 54454.000
 Mean: 55090.037
 Std. dev.: 16554.377
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
191: 27318.0
187: 30573.0
202: 36960.0
207: 37233.0
197: 39628.0
177: 40244.0
147: 42508.0
172: 43010.0
198: 44468.0
201: 47403.0
206: 49154.0
186: 50631.0
188: 53428.0
203: 54454.0
196: 56691.0
178: 57861.0
208: 59948.0
192: 60098.0
183: 60631.0
173: 61142.0
182: 62983.0
193: 63911.0
146: 65130.0
176: 77668.0
171: 78436.0
148: 86061.0
181: 99859.0

USE THESE NUMBERS TO SET COREDIVE ALPHA RAREFACTION LIMIT (MEAN - 2*SD) _ arround 22k in this case... 
could guess from multiqc read counts and typical bacterial or fungal flash rates

EIGENVALUES WARNINGS
PYTHON ELEMENT-WISE VS SCALAR COMPARISON

/home/mark/miniconda3/envs/qiime1/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.00688267504561 and the largest is 0.620990940226.
  RuntimeWarning
/home/mark/miniconda3/envs/qiime1/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
/home/mark/miniconda3/envs/qiime1/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
/home/mark/miniconda3/envs/qiime1/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
  
  
ITS1F (5’- CTT GGT CAT TTA GAG GAA GTA A-3’) and 58A2R (5’-CTG CGT TCT TCA TCG AT-3’) for amplification of fungi

what to do about NA and uncertain rank asignemnts?

empirical rules for picking top N and k/a?
what are the rnak and k/a for the samples that shwo a differenrce?

191
pre flash 29.5k
post flash 29.1k
post split & pick 27260.0 (93.8%)

problem:
/home/mark/gitstage/uderica/nordqiime/15_Hil_CS/Fungi/ITS/split/178/split_library_log.txt
/home/mark/gitstage/uderica/nordqiime/15_Hil_CS/Fungi/ITS/extendedFrags/178.extendedFrags.fastq



http://bioinformatics.cvr.ac.uk/blog/short-command-lines-for-manipulation-fastq-and-fasta-sequence-files/

/home/mark/Downloads/ncbi-blast-2.6.0+/bin/makeblastdb -in ~/Downloads/yeast.nt -parse_seqids -dbtype nucl

find a control that is internally consistent and a treated sample that is consistent by\uyt distanc\t from the cntrol

https://github.com/weizhongli/cdhit

# fastq to fasta... I have been happy with https://www.biostars.org/p/85929/
http://seqanswers.com/forums/showthread.php?t=6888
https://www.biostars.org/p/85929/
sed -n '1~4s/^@/>/p;2~4p' 191.notCombined_1.fastq > 191.notCombined_1.fasta

https://github.com/marbl/Krona/wiki/Importing-from-supported-tools
https://www.biostars.org/p/52623/

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
https://www.ncbi.nlm.nih.gov/books/NBK279675/

# discussion about mega for identifying origi of contaminants
https://www.biostars.org/p/54081/

# qiime script that takes R1 and R2 reads (wouldn't require FLASH)
http://qiime.org/tutorials/rtax.html


http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# DESeq2 multie factors
https://support.bioconductor.org/p/66263/


http://biocore.github.io/pynast/

# Aquatic environmental DNA detects seasonal fish abundance and habitat preference in an urban estuary
# Mark Y. Stoeckle , Lyubov Soboleva, Zachary Charlop-Powers
http://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0175186

# contamination discussion (PCR improvements)
https://www.researchgate.net/post/Pyrosequencing_of_a_bacterial_community_from_plant_samples

http://qiime.org/scripts/filter_taxa_from_otu_table.html
https://rdrr.io/bioc/phyloseq/man/prune_taxa-methods.html

https://www.researchgate.net/post/How_long_is_the_ITS_internal_transcribed_spacer_PCR_fragment_of_Tulasnella_fungi


looking for contaminants

cdhit... how to get counts/cluster?

having trouble compiling/using:
cdhit otu
cdhit dup

usearch?

add alpha to diffy abundance plot

mito/chloro


