# STAR2count_by_starcode.py


### Usage Example

`python /home/yuanhao/single_cell/scripts/STAR2count_by_starcode.py --gene-tag=GE \\
  --outsam "/home/yuanhao/single_cell/GSE110823_output/SRR6750041/star_results/corrected.sam" \\
  --outcount "/home/yuanhao/single_cell/GSE110823_output/SRR6750041/star_results/count.tsv" \\
  --starcode-path "/home/yuanhao/softwares/starcode/starcode" \\
  "/home/yuanhao/single_cell/GSE110823_output/SRR6750041/star_results/test.sam"`






As the article "Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding"

http://science.sciencemag.org/content/early/2018/03/14/science.aam8999 

says:
  "We then used Starcode (45) to collapse UMIs of aligned reads that were within 1 nt mismatch of another UMI, assuming the two aligned reads were also from the same UBC. Each original barcoded cDNA molecule is amplified before tagmentation and subsequent PCR, so a single UMI-UBC combination can have several distinct cDNA reads corresponding to different parts of the transcript. Occasionally STAR will map these different reads to different genes. As a result, we chose the most frequently assigned gene as the mapping for the given UMI-UBC combination. We then generated a matrix of gene counts for each cell (N x K matrix, with N cells and K genes)."
 
I tried to follow it and write this script from starcode-umi to facilitate this procedure
 
  - starcode is from https://github.com/gui11aume/starcode
 
  - STAR is from https://github.com/alexdobin/STAR
 
  - dropseq_tools is from http://mccarrolllab.com/dropseq/
 
  - umi_tools is from https://github.com/CGATOxford/UMI-tools
 
### Input

The inout file is a sam file from TagReadWithGeneExon, which is extracted UMI and barcodes by umi_tools and aligned by STAR sequentially.
 
The input samfile format maybe:
 
`SRR6750041.55844198_AACCACATACACAGAATGGAACAACAGATC_ATCTCGTGTG   0       GRCh38_19       8455426 3       66M     *       0       0       GAATGAGAAGAGGAAGGAGAAAAACATAAAAAGAGGAGGCAATCGCTTTGAGCCATATGCCAACCC    AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEA      GE:Z:GRCh38_HNRNPM   XF:Z:CODING      NH:i:2  HI:i:1  nM:i:1  AS:i:62 GS:Z:+
SRR6750041.54895966_ACTATGCAACACAGAATGGAACAACAGATC_CGTCATGTAT   0       GRCh38_19       8455426 3       66M     *       0       0       GAATGAGAAGAGGAAGGAGAAAAACATAAAAAGAGGAGGCAATCGCTTTGAGCCATATGCCAACCC    AAAAAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEAEEEE66      GE:Z:GRCh38_HNRNPM   XF:Z:CODING      NH:i:2  HI:i:1  nM:i:1  AS:i:62 GS:Z:+
SRR6750041.55054756_AACCACATACACAGAATGGAACAACAGATC_ATCTCGTGTG   0       GRCh38_19       8455426 3       66M     *       0       0       GAATGAGAAGAGGAAGGAGAAAAACATAAAAAGATGAGGCAATCGCTTTGAGCCATATGCCAACCC    AAAAAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEAEAEEEAEAEEE      GE:Z:GRCh38_HNRNPM   XF:Z:CODING      NH:i:2  HI:i:1  nM:i:2  AS:i:60 GS:Z:+`


  - In which AACCACATACACAGAATGGAACAACAGATC is cell barcode and ATCTCGTGTG is UMI
 
### Reference:
 
  1.Rosenberg A B, Roco C M, Muscat R A, et al. Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding[J]. Science, 2018, 360(6385): 176-182.
 
  2.Zorita E, Cusco P, Filion G J. Starcode: sequence clustering based on all-pairs search[J]. Bioinformatics, 2015, 31(12): 1913-1919.
 
  3.Macosko E Z, Basu A, Satija R, et al. Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets[J]. Cell, 2015, 161(5): 1202-1214.
 
  4.Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21.

 
 
