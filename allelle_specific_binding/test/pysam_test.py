import pysam
bamfile="/mnt/isilon/sfgi/dudekm/raw_data/brandon_liver_ATAC/chr_renamed_bams/3782-BW-99.bam"
samfile = pysam.AlignmentFile(bamfile, "rb")
chrom = "chr10"
pos = 718264

unpaired_outfilename = "bams/unpaired_test.bam"
paired_outfilename = "bams/paired_test.bam"

unpaired_outfile = pysam.AlignmentFile(unpaired_outfilename, "wb", template=samfile)
paired_outfile = pysam.AlignmentFile(paired_outfilename, "wb", template=samfile)

for pileup in samfile.pileup(chrom, pos-1, pos):
    if pileup.reference_pos == pos-1:
            for read in pileup.pileups:
                    print(read.alignment.query_name)
                    print("is.read1 =", read.alignment.is_read1)
                    print("is.read2 =", read.alignment.is_read2)
                    mate = samfile.mate(read.alignment)
                    print(mate.query_name)
                    print("mate is.read1 =", mate.is_read1)
                    print("mate is.read2 =", mate.is_read2)
                    unpaired_outfile.write(read.alignment)
                    paired_outfile.write(read.alignment)
                    paired_outfile.write(mate)

unpaired_outfile.close()
paired_outfile.close()

pysam.sort("-o", unpaired_outfilename + ".sorted", unpaired_outfilename)
pysam.sort("-o", paired_outfilename + ".sorted", paired_outfilename)

unpaired_outfilename += ".sorted"
paired_outfilename += ".sorted"

pysam.index(unpaired_outfilename)
pysam.index(paired_outfilename)


# Check test outfiles
print("\n---\nChecking test outfiles\n---")
unpaired_outfile = pysam.AlignmentFile(unpaired_outfilename, "rb")
paired_outfile = pysam.AlignmentFile(paired_outfilename, "rb")

print("Unpaired...")
for pileup in unpaired_outfile.pileup(chrom, pos-1, pos):
    if pileup.reference_pos == pos-1:
            for read in pileup.pileups:
                    print(read.alignment.query_name)
                    print("is.read1 =", read.alignment.is_read1)
                    print("is.read2 =", read.alignment.is_read2)

print("paired...")
for pileup in paired_outfile.pileup(chrom, pos-1, pos):
    if pileup.reference_pos == pos-1:
            for read in pileup.pileups:
                    print(read.alignment.query_name)
                    print("is.read1 =", read.alignment.is_read1)
                    print("is.read2 =", read.alignment.is_read2)
                    
