from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd

for record in SeqIO.parse("../syn3_9_01.fa", "fasta"):
    seq = str(record.seq)
    dic = pd.read_csv("../PCRtag_yeast_chr03_3_40.csv", header=0)
    my_seq = Seq(seq, IUPAC.unambiguous_dna)
    print(seq)
    print(my_seq.complement())
    print(my_seq.reverse_complement())
