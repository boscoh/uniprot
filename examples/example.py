#!/usr/bin/env python
"""Example usage of the uniprot module."""

import os
import pprint
import uniprot

EXAMPLES_DIR = os.path.dirname(os.path.abspath(__file__))
TESTS_DATA_DIR = os.path.join(os.path.dirname(EXAMPLES_DIR), 'tests', 'data')

os.system('rm -f cache* 2>/dev/null')

# Example 1 - reading a fasta file
fasta_file = os.path.join(EXAMPLES_DIR, 'example.fasta')
seqids, fastas = uniprot.read_fasta(fasta_file)
print("Example 1 - Reading FASTA file:")
pprint.pprint(seqids, indent=2)
print()

# Example 2 - map identifiers for RefSeq to Uniprot
seqids = "NP_000508.1  NP_001018081.3".split()
pairs = uniprot.batch_uniprot_id_mapping_pairs(
    'RefSeq_Protein', 'UniProtKB', seqids)
print("Example 2 - RefSeq to UniProt mapping:")
pprint.pprint(pairs, indent=2)
print()

# Example 3 - get UniProt metadata
uniprot_seqids = [j for i, j in pairs]
uniprot_data = uniprot.batch_uniprot_metadata(uniprot_seqids, 'cache')
print("Example 3 - UniProt metadata:")
pprint.pprint(uniprot_data, indent=2)
print()

# Example 4 - parse for isoforms in metadata
if os.path.exists('cache/metadata.0.txt'):
    text = open('cache/metadata.0.txt').read()
    uniprot_data = uniprot.parse_isoforms(text)
    print("Example 4 - Isoform parsing:")
    pprint.pprint(uniprot_data)
    print()

# Example 5 - chaining commands to map seqids
seqids = "EFG_MYCA1 YP_885981.1 ENSG00000196176 Q91ZU6-8".split()
uniprot_data = uniprot.get_metadata_with_some_seqid_conversions(seqids, 'cache2')
print("Example 5 - Metadata with seqid conversions:")
pprint.pprint(uniprot_data, indent=2)
print()

# Example 6 - verify isoform sequence
isoform_fasta = os.path.join(TESTS_DATA_DIR, 'isoform', 'Q91ZU6-8.fasta')
if os.path.exists(isoform_fasta):
    read_seqids, fastas = uniprot.read_fasta(isoform_fasta)
    test_sequence = list(fastas.values())[0]['sequence']
    if 'Q91ZU6-8' in uniprot_data:
        print("Example 6 - Isoform sequence match:", 
              test_sequence == uniprot_data['Q91ZU6-8']['sequence'])
        print()

# Example 7 - check isoforms
isoform_txt = os.path.join(TESTS_DATA_DIR, 'isoform', 'Q91ZU6.txt')
if os.path.exists(isoform_txt):
    seqids = ["Q91ZU6-{}".format(i) for i in [1, 2, 3, 4, 5, 6, 8]]
    txt = open(isoform_txt).read()
    results = uniprot.parse_uniprot_metadata_with_seqids(seqids, txt)
    print("Example 7 - Isoform sequence verification:")
    for seqid in seqids:
        fasta_file = os.path.join(TESTS_DATA_DIR, 'isoform', seqid + '.fasta')
        if os.path.exists(fasta_file):
            read_seqids, fastas = uniprot.read_fasta(fasta_file)
            test_sequence = list(fastas.values())[0]['sequence']
            print(seqid, test_sequence == results[seqid]['sequence'])

