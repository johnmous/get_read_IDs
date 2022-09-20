import pysam
import multiprocessing
import argparse
import sys
import time

parser = argparse.ArgumentParser(description="Get list of read IDs from alignment file")
parser.add_argument("--aln",
                    dest="aln_file",
                    help="Alignment file")
parser.add_argument("--brc",
                    dest="brc_list",
                    help="File with barcodes, one barcode per line")
parser.add_argument("--cores",
                    dest="cores",
                    type=int,
                    help="Number of cores available")
parser.add_argument("--out_file",
                    dest="out_file",
                    help="File to save the read IDs")
args = parser.parse_args()

# Open alignment file
samfile = pysam.AlignmentFile(args.aln_file, "rb")
chroms = samfile.references
samfile.close()


# Load the cell brcode list
with open(args.brc_list) as file:
    cell_bc = file.read().splitlines()
cell_bc.pop(0)


def fetch_ids(chrom):
    """
    Loop through records and select those with corrected barcode tag (CB) present
    in the list of provided barcodes
    :param chrom: name of chromosome to work on
    :return: list of read IDs
    """
    start = time.time()
    read_ids = []
    samfile = pysam.AlignmentFile(args.aln_file, "rb")
    for read in samfile.fetch(chrom):
        try:
            if read.get_tag("CB") in cell_bc:
                read_ids.append(read.query_name)
        except KeyError:
            pass
    samfile.close()
    end = time.time()
    print(f"Finished processing chromosome: {chrom}. It took {end - start} seconds to calculate.", file=sys.stderr)
    return read_ids


if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=args.cores)
    outputs = pool.map(fetch_ids, chroms)

    # outputs is a list of lists, flatten it
    ids = [item for sublist in outputs for item in sublist]

    # Convert to set to remove possible duplicates, then save to file
    read_ids = set(ids)
    with open(args.out_file, 'w') as f:
        for i in read_ids:
            f.write(f"{i}\n")
