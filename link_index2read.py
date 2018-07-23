#!/usr/bin/env python

import argparse
import sys
import gzip
import time

# Parse arguments
parser = argparse.ArgumentParser(description='Link index to reads')
parser.add_argument('--out-file', required=True, type=str, help="output directory")
parser.add_argument('--out-file2', required=False, default="", type=str, help="output directory")
parser.add_argument('--index', required=True, type=str, help="input index file")
parser.add_argument('--index-base_number', required=False, default=6, type=int, help="input index file")
parser.add_argument('--read-file2', required=False, default="", type=str, help="input read file 2")

parser.add_argument('read_file', help="input read fastq")


params = parser.parse_args()


def make_read_form_output(index_line, fastq_form, read_f):
    read_form_output = ""
    for line_unit in range(fastq_form):
        read_line = read_f.readline()
        if line_unit == 0:
            read_line_part0 = read_line.rsplit(":", 1)[0]
            read_line = "{}:{}\n".format(read_line_part0, index_line[:params.index_base_number])
        read_form_output = "{}{}".format(read_form_output, read_line)
    return read_form_output

if __name__ == '__main__':
    line_number = 0
    fastq_form = 4



    if params.read_file2 != "":
        with gzip.open(params.index) as index_f,\
             gzip.open(params.read_file) as read1_f,\
             gzip.open(params.read_file2) as read2_f, \
             open(params.out_file, "w") as out1_f, \
             open(params.out_file2, "w") as out2_f:
            read_index = 0
            while not index_f.readline() == "":
                for line_unit in range(fastq_form - 1):
                    if line_unit == 0:
                        index_line = index_f.readline()
                    else:
                        _ = index_f.readline()
                read_index += 1
                read1_form_output = make_read_form_output(index_line, fastq_form, read1_f)
                read2_form_output = make_read_form_output(index_line, fastq_form, read2_f)
                out1_f.write(read1_form_output)
                out2_f.write(read2_form_output)
                if read_index % 100000 == 0:
                    print("process {} reads".format(read_index))

    else:
        with gzip.open(params.index) as index_f,\
             gzip.open(params.read_file) as read1_f,\
             open(params.out_file, "w") as out1_f:
            read_index = 0
            while not index_f.readline() == "":
                for line_unit in range(fastq_form - 1):
                    if line_unit == 0:
                        index_line = index_f.readline()
                    else:
                        _ = index_f.readline()
                read_index += 1
                read1_form_output = make_read_form_output(index_line, fastq_form, read1_f)
                out1_f.write(read1_form_output)
                if read_index % 100000 == 0:
                    print("process {} reads".format(read_index))
