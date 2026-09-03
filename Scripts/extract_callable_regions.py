from Bio import SeqIO

input_fasta = "CocciRef_GCA_000149335.2.masked.fna"
output_bed = "callable_regions.bed"

with open(output_bed, "w") as bed_out:
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper()
        start = None
        for i, base in enumerate(seq):
            if base != "N":
                if start is None:
                    start = i
            else:
                if start is not None:
                    bed_out.write(f"{record.id}\t{start}\t{i}\n")
                    start = None
        # Close any open region at end of sequence
        if start is not None:
            bed_out.write(f"{record.id}\t{start}\t{len(seq)}\n")
