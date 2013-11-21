import os
import uniprot
import pprint


# Clean up caches
os.system('rm cache*')

# Example 1 - reading a fasta file
seqids, fastas = uniprot.read_fasta('example.fasta')
pprint.pprint(seqids, indent=2)

# Example 2 - map identifiers for RefSeq to Uniprot
seqids = "NP_000508.1  NP_001018081.3".split()
pairs = uniprot.batch_uniprot_id_mapping_pairs(
  'P_REFSEQ_AC', 'ACC', seqids)
pprint.pprint(pairs, indent=2)

# Example 2 - get UniProt metadata
uniprot_seqids = [j for i,j in pairs]
uniprot_data = uniprot.batch_uniprot_metadata(
    uniprot_seqids, 'cache')
pprint.pprint(uniprot_data, indent=2)

# Example 3 - parse for isoforms in metadata
text = open('cache.0.txt').read()
uniprot_data = uniprot.parse_isoforms(text)
pprint.pprint(uniprot_data)

# Example 4 - chaining commands to map seqids
seqids = "EFG_MYCA1 YP_885981.1 ENSG00000196176".split()
uniprot_data = uniprot.get_metadata_with_some_seqid_conversions(
    seqids, 'cache2')
pprint.pprint(uniprot_data, indent=2)

