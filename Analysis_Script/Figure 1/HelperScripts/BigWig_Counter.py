import pyBigWig
from timeit import default_timer as clock
from multiprocess import Pool
import pandas as pd


def fetch_counts(args):
    """
    Gets a vector of counts for the bed_regions for one individual bw_file, so that we can parallelize over the
    bw files in the outer function.
    """
    bw_file, bed_regions = args
    collected_counts = []
    collected_errors = []
    this_bw = pyBigWig.open(bw_file)

    # Check if the bigwig uses the 'chr' prefix or not. The bed_regions always have the chr-prefix added.
    this_bw_chroms = this_bw.chroms()
    if next(iter(this_bw_chroms.keys())).startswith('chr'):
        remove_prefix = ''
    else:  # If the bigwig doesn't have the prefix we remove it again.
        remove_prefix = 'chr'

    for region in bed_regions:
        try:
            bw_count = this_bw.stats(region[0].replace(remove_prefix, ''), int(region[1]), int(region[2]), type="mean")[0]
            if bw_count is None:
                bw_count = 0
                collected_errors.append([region, "count was NaN"])
            collected_counts.append(bw_count)
        except (RuntimeError, IndexError):
            collected_counts.append(0)
            collected_errors.append([region, "Invalid interval bounds"])

    if not len(collected_counts) == len(bed_regions):
        print("ERROR: Mismatch between fetched counts and number of bed regions", bw_file)
    return collected_counts, collected_errors


def bigwig_counts(bed_file, bigwigs, n_cores=1):
    """
    For a bed file or BedTool object gets the mean signal for each bigwig file. If pyBigWig can't retrieve a count it
    will be set to 0, and the region listed in the error list.

    Args:
        bed_file: Either a BedTools object or a path to a bed-file.
        bigwigs: List/dict of bigwigs for which pyBigWig will be used to extract the mean signal for. If a list of files is given, the column in the resulting DataFrame will have the full file path. If a dict {key: file} is given, the columns will be the keys.

    Returns:
        tuple:
            - **region_counts**: Will give a dataframe of the bed-file's regions with a column for each bigwig file and the file name as column.
            - **errors**: Dictionary with an entry for each biwgig file holding a DataFrame with the regions that are also part of the bed_regions but failed due to not returning a count or throwing an error of invalid interval bounds. Those have a count of 0. Note, the errors are rather unreliable, the behaviour of the pyBigWig is a bit elusive.
    """
    start = clock()

    if type(bed_file) == str:
        bed_regions = [x.strip().split('\t') for x in open(bed_file).readlines() if not x.startswith('#')]
    else:
        bed_regions = [str(x).strip().split('\t') for x in str(bed_file).strip().split('\n')]

    bigwig_order = bigwigs
    df_columns = bigwig_order
    if type(bigwigs) == set:
        bigwig_order = df_columns = list(bigwigs)
    if type(bigwigs) == dict:
        bigwig_order = list(bigwigs.values())
        df_columns = list(bigwigs.keys())

    process_pool = Pool(processes=n_cores)
    bw_counts = process_pool.map(fetch_counts, map(lambda x: (x, bed_regions), bigwig_order))
    process_pool.close()
    region_counts = pd.DataFrame([x[0] for x in bw_counts]).T
    region_counts = pd.concat([pd.DataFrame([x[:3] for x in bed_regions]), region_counts], ignore_index=True, axis=1)
    region_counts.columns = ['#chr', 'start', 'end'] + df_columns

    errors = {}
    for i_bw, bw_errors in enumerate([x[1] for x in bw_counts]):
        if bw_errors:
            errors[df_columns[i_bw]] = pd.DataFrame(bw_errors, columns=['region', 'error'])
    if errors:
        print("WARNING: some regions returned errors, see returned errors df.")

    print('Bigwig counts fetched', clock() - start)
    return region_counts, errors


