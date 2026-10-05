from collections import Counter
import pandas as pd
import numpy as np
import random
from timeit import default_timer as clock
import pybedtools
import sys
from pybedtools import BedTool
import gzip
from multiprocess import Pool
import copy
import scipy.stats
import statsmodels.stats.multitest
from sklearn.decomposition import PCA
import os
import fnmatch
import re
from itertools import chain
import GTF_Processing

"""Collection of functions that didn't fit well anywhere else. Also includes functions which are
required by multiple other modules."""


def sanitize_path(path_string):
    """Function to remove unwanted characters from a file path before saving. Put here since we need it for a lot
    of functions and eases adjustment e.g. for OS."""
    return re.sub(r'[^A-Za-z0-9._\-\\\\////]+', '', path_string)
    # return path_string.replace(' ', '').replace(':', '').replace("°", "Degree").replace()


def open_genes_multibed(bed_dict, annotation):
    """For each of the bed files in bed_dict {tag: bed file} returns the set of genes that intersect with at least
    one of their promoters with a peak."""
    promoter_bed = GTF_Processing.gene_window_bed(annotation, extend=200, tss_type='all')
    open_genes = {}
    for tag, bed in bed_dict.items():
        print(bed)
        if type(bed) == str:
            bed = BedTool(bed)
        open_genes[tag] = set([x.fields[3] for x in promoter_bed.intersect(bed, u=True)])

    return open_genes


def open_genes(activity_file, annotation, first_col, tmp_dir):
    """
    NOTE: The GTF_Processing.gene_window_bed fulfils this function when not working with an activity column but just a bed file.
    Takes a bed-styled file with the activities of regions across multiple samples/conditions in columns.
    Takes all annotated promoters ±200bp for genes and checks which genes have an intersecting peak that has an
    activity > 0. The activity file requires a header preceded by #.
    @param open_genes: A dictionary where the keys are the column headers from the activity file and the
    values the set of genes with a nonzero peak at one of their promoters in that column.
    """
    pybedtools.set_tempdir(tmp_dir)
    start_prom = clock()
    promoter_bed = GTF_Processing.gene_window_bed(annotation, extend=200, tss_type='all', open_regions=activity_file)
    print(clock() - start_prom, "Got promoter bed")
    open_genes = {c: set() for c in open(activity_file).readline().strip().split('\t')[first_col:]}
    # Get the index of the columns in the bed intersection. The promoter columns are the first 4 entries.
    col_idx = {i: c for i, c in enumerate(open(activity_file).readline().strip().split('\t')[first_col:])}

    helper = 0
    for chro in ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']:
        print(chro)
        peak_nonzeros = {}
        start_chro = clock()
        bed_chro = ''
        with open(activity_file, 'r') as bed_in:
            bed_in.readline()  # Skip header.
            for entry in bed_in:
                if entry[:entry.find('\t')] == chro:
                    entry = entry.split('\t')
                    peak_nonzeros['\t'.join(entry[1:3])] = [c for c, val in enumerate(entry[first_col:]) if
                                                            float(val) > 0]
                    bed_chro += '\t'.join(entry[:3]) + '\n'
        print(clock() - start_chro, "bed fetched")
        chr_promoter = BedTool(''.join([str(x) for x in promoter_bed if x.fields[0] == chro]), from_string=True)
        for inter in chr_promoter.intersect(BedTool(bed_chro, from_string=True), wo=True):
            helper += 1
            if helper % 10000 == 0:
                print(helper)
            gene = inter.fields[3]
            for c in peak_nonzeros['\t'.join(inter.fields[5:7])]:
                open_genes[col_idx[c]].add(gene)
        print(clock() - start_chro)


def fn_patternmatch(pattern):
    """
    Grabs all files in the file system that match the pattern and returns a dictionary with {file: wildcard}.
    Only works if the wildcard is in the file name and not in a directory name.
    Files that start with ._ (e.g. temporary or system files) will be skipped.
    E.g. for a directory BirdCollection/ that contains Bird_Kakapo.txt and Bird_Kea.txt:
    fn_patternmatch("BirdCollection/*.txt" = {"BirdCollection/Bird_Kakapo.txt": "Kakao", "BirdCollection/Bird_Kea.txt": "Kea"}
     """
    parent_folder = '/'.join(pattern.split('/')[:-1]) + '/'
    children_pattern = pattern.split('/')[-1]
    re_pattern = re.compile(children_pattern.replace('*', '(.*)'))
    matched_files = {parent_folder + x: re_pattern.search(x).group(1)
                     for x in os.listdir(parent_folder)
                     if fnmatch.fnmatch(x, children_pattern) and not x.startswith('._')}
    return matched_files
