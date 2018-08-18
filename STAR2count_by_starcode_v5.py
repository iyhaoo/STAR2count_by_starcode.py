#!/usr/bin/env python

import argparse
import sys
import subprocess
from collections import Counter
import time
import copy

# Parse arguments
parser = argparse.ArgumentParser(description='Use starcode to demultiplex reads containing UMI identifiers.')
# parser.add_argument('--umi-len', required=True, type=int, help="UMI lenght (must be at the beginning of the read)")
parser.add_argument('--out-count', required=False, default="", type=str, help="output expression matrix")
parser.add_argument('--out-sam', required=False, default="", type=str, help="tem_sam")
parser.add_argument('--samtools-path', required=False, default="", type=str, help="use samtools to read bam-file")
parser.add_argument('--umi-d', required=False, default=0, type=int, help="UMI match distance (default: 0)")
parser.add_argument('--gene-tag', required=True, type=str, help="sam gene tag")
parser.add_argument('--umi-cluster', required=False, default='mp', type=str,
                    help="Algorithm used to cluster UMIs: 'mp' for message passing, 's' for spheres or 'cc' for connected components (default: 'mp')")
parser.add_argument('--umi-threads', required=False, default=1, type=int,
                    help="Starcode threads for UMI clustering (recommended: 1, default: 1)")
parser.add_argument('--umi-cluster-ratio', required=False, default=3, type=int,
                    help="Min size ratio for merging clusters in UMI message passing (default: 3)")
parser.add_argument('--starcode-path', required=False, default='./starcode', type=str,
                    help="Path to starcode binary file (default: './starcode')")
parser.add_argument('sam_file',
                    help="Input file in sam format.")


params = parser.parse_args()

# Starcode options
path = params.starcode_path
umi_cluster = params.umi_cluster
umi_ratio = params.umi_cluster_ratio
samtool_path = params.samtools_path

# UMI and SEQ params
umi_tau = params.umi_d  # umi mismatches

#count_options
gene_tag = params.gene_tag


if __name__ == '__main__':
    # Check arguments
    if params.out_count != "":
        with open(params.out_count, "w") as f:
            f.write("")
    clust_opts = ['s', 'mp', 'cc']
    if umi_cluster not in clust_opts:
        sys.stderr.write("Error, cluster option must be 's', 'mp' or 'cc'.\n")
        sys.exit(1)
    # UMI options
    if umi_cluster == 's':
        umi_args = "-s"
    elif umi_cluster == 'mp':
        umi_args = "-r{}".format(umi_ratio)
    elif umi_cluster == 'cc':
        umi_args = "-c"
    umiarg = [path, "--seq-id", "-qd" + str(umi_tau), umi_args, "-t" + str(params.umi_threads)]
    # Subprocess: Cluster UMI.
    umiproc = subprocess.Popen(umiarg, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

    sys.stderr.write("parsing barcodes\n")
    lastTime = time.time()
    tmp_file = params.sam_file.replace("\\", "/").rsplit("/", 1)[0] + "/s2ctmp"

    if params.sam_file.rsplit(".", 1)[1] == "bam":
        is_bam = True
        if samtool_path == "":
            sys.stderr.write("Error, cannot read bam-file without samtools.\n")
            sys.exit(1)
        samtools_arg = [samtool_path, "view", "-h", params.sam_file]
        samtools_proc = subprocess.Popen(samtools_arg, stdout=subprocess.PIPE)
        f = samtools_proc.stdout
        f_c = copy.deepcopy(f)
    else:
        is_bam = False
        f = open(params.sam_file)
        f_c = open(params.sam_file)
    for l_i, line in enumerate(f):
        if line[0] != "@":
            # Starcode seqs and UMIs independently.
            line = line.rstrip()
            split_line = line.split("\t")
            read_id_region = split_line[0]
            split_read_id_region = read_id_region.split("_", 2)
            umi_s = split_read_id_region[2]
            umiproc.stdin.write(umi_s + "\n")
        if l_i % 1000000 == 0:
            sys.stderr.write("Read {} lines in {} seconds\n".format(l_i, time.time() - lastTime))
            lastTime = time.time()
        # if l_i > 300000:
        #    break ###test
    if not is_bam:
        f.close()
    # Close pipes and let starcode run.
    umiproc.stdin.close()
    sys.stderr.write("Start make dictionary (use {} seconds)\n".format(time.time() - lastTime))
    lastTime = time.time()
    mapping_dict = {}
    for line in umiproc.stdout:  # from 1 to len()
        line_split = line.split("\t", 2)
        for index in line_split[2].replace("\n", "").split(","):
            mapping_dict[index] = line_split[0]
    sys.stderr.write("Find {} reads\n".format(len(mapping_dict.keys())))
    if params.out_sam != "":
        write_collapsed_sam = True
        o_f = open(params.out_sam, "w")
    else:
        write_collapsed_sam = False

    barcode_dict = {}
    for l_i, line in enumerate(f_c):
        if line[0] != "@":
            into_star_code_line = l_i + 1 # from 1 to len()
            line = line.rstrip()
            split_line = line.split("\t")
            read_id_region = split_line[0]
            split_read_id_region = read_id_region.split("_", 2)
            barcode_s = split_read_id_region[1]
            tags = [x.split(":", 1)[0] for x in split_line[11:]]  # ##after quality
            if params.gene_tag in tags:
                out_str = "\t".join(["_".join(split_read_id_region[:2] + [mapping_dict[str(into_star_code_line)]])] + split_line[1:])
                try:
                    try:
                        barcode_dict[barcode_s][mapping_dict[str(into_star_code_line)]] = "{}\t{}".format(
                            barcode_dict[barcode_s][mapping_dict[str(into_star_code_line)]],
                            split_line[11 + tags.index(gene_tag)].rsplit(":", 1)[1])
                    except:
                        barcode_dict[barcode_s][mapping_dict[str(into_star_code_line)]] = \
                            split_line[11 + tags.index(gene_tag)].rsplit(":", 1)[1]
                except:
                    barcode_dict[barcode_s] = {}
                    barcode_dict[barcode_s][mapping_dict[str(into_star_code_line)]] = \
                        split_line[11 + tags.index(gene_tag)].rsplit(":", 1)[1]
                if l_i % 1000000 == 0:
                    sys.stderr.write("Processing {} reads with gene_tag in {} seconds\n".format(l_i, time.time() - lastTime))
                    lastTime = time.time()
                if write_collapsed_sam:
                    o_f.write("{}\n".format(out_str))
    if not is_bam:
        f_c.close()
    if write_collapsed_sam:
        o_f.close()

    out_dict = {}
    for cell, umis in barcode_dict.items():
        for gene_express in umis.values():
            genelist = gene_express.split("\t")
            # we chose the most frequently assigned gene as the mapping for the given UMI-UBC combination
            canonical_gene = Counter(genelist).keys()[0]
            gene_cell_element = "{}\t{}".format(canonical_gene, cell)
            try:
                out_dict[gene_cell_element] += 1
            except:
                out_dict[gene_cell_element] = 1
    sys.stderr.write("Got {} unique UMIs\n".format(sum(out_dict.values())))
    if params.out_count != "":
        with open(params.out_count, "w") as f:
            for key, val in out_dict.items():
                f.write("{}\t{}\n".format(key, val))
    else:
        for key, val in out_dict.items():
            print("{}\t{}".format(key, val))