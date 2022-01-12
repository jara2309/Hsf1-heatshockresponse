import sys

chr_dict = {"chrX": "NC_003284.9", "chrM": "NC_00328.1", "chrV": "NC_003283.11", "chrIV": "NC_003282.8",
            "chrIII": "NC_003281.10", "chrII": "NC_003280.10",
            "chrI": "NC_003279.8"}  # dictionary to rename the .bed12 file

downloaded_bed = sys.argv[1]
new_bed = sys.argv[2]

new_bed = open(new_bed, "w+")

with open(downloaded_bed, "r") as f:
    for line in f:
        split_line = line.rstrip().split()
        if split_line:
            chr_name = split_line[0]  # which column to chose    bed 3

            new_bed.write(line.replace(chr_name, chr_dict.get(chr_name)))  # replace that with the new name
new_bed.close()
