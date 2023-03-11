import csv


def write_tsv(input, handle):
    tsv_writer = csv.writer(
        handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    for i in input:
        tsv_writer.writerow(i)
