
from pprint import pprint
import re
import uniprot





text = open('Q91ZU6.txt').read()
uniprot_data = uniprot.parse_isoforms(text)
for entry in uniprot_data.values():
  isoforms = entry['isoforms']
  var_seqs = entry['var_seqs']
  for isoform in isoforms.values():
    fasta_db = isoform['seqid'] + '.fasta'
    seqids, fastas = uniprot.read_fasta(fasta_db)
    fasta_sequence = fastas.values()[0]['sequence']
    assert fasta_sequence == isoform['sequence']
