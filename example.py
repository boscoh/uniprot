import os
import uniprot
import pprint


# Clean up caches

os.system('rm *output* *cache*')


# Example 1 - reading a fasta file

seqids, fastas = uniprot.read_fasta('example.fasta')
pprint.pprint(seqids, indent=2)


# Example 2 - batch read identifier mappings with
# prespecified identifier types

seqids = """
NP_000508.1  NP_001018081.3
""".split()

pairs = uniprot.batch_uniprot_id_mapping_pairs(
  'P_REFSEQ_AC', 'ACC', seqids)

pprint.pprint(pairs, indent=2)


# Example 3 - sequential identifier mapping to UniProt 
# identifiers using robust but slow method

seqids = """
EFG_MYCA1 YP_885981.1 CpC231_1796
""".split()

mapping = uniprot.sequentially_convert_to_uniprot_id(
    seqids, 'cache.json')

uniprot_seqids = mapping.values()


# Example 4 - get UniProt metadata

uniprot_data = uniprot.batch_uniprot_metadata(
    uniprot_seqids, 'cache2.txt')

pprint.pprint(uniprot_data, indent=2)

for l in open('cache2.txt'):
  print l.strip()

uniprot.write_fasta('example.output.fasta', uniprot_data, uniprot_seqids)


# Example 5 - chaining commands to make your own
# special mapper


def map_to_refseq(seqids):
  uniprot_mapping = uniprot.sequentially_convert_to_uniprot_id(
      seqids, 'func.cache.json')
  uniprot_ids = uniprot_mapping.values()
  pairs = uniprot.batch_uniprot_id_mapping_pairs(
    'ACC', 'P_REFSEQ_AC', uniprot_ids)
  mapping = {}  
  for seqid in seqids:
    if seqid in uniprot_mapping:
      uniprot_id = uniprot_mapping[seqid]
    for pair in pairs:
      if uniprot_id == pair[0]: 
        mapping[seqid] = pair[1]
  os.remove('func.cache.json')
  return mapping


seqids = """
EFG_MYCA1 YP_885981.1 CpC231_1796
""".split()

mapping = map_to_refseq(seqids)

pprint.pprint(mapping, indent=2)