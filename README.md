

# UNIPROT.PY


`uniprot.py` provides a Python interface to the UniProt website <http://uniprot.org> to: 

1. Map between types of protein sequence identifiers

2. Fetch metadata for proteins such as organism, the 
sequence and GO annotations.


## Installation

`uniprot.py` has one dependency: the brilliant `requests.py` that wraps the dire `httplib` in the standard Python library. To install:

    >> pip install uniprot


## Examples

Before you run any of the examples, import the module:

    import uniprot

I also like importing pprint to print out dictionaries,
in order to see what I am getting back:

    import pprint

A convenience function is provided to read seqids and sequences from a fasta file:

    seqids, fastas = uniprot.read_fasta('example.fasta')


## Fetching identifier mappings

If we have a bunch of sequence identifiers that belong to one of the known 
seqid types understood by Uniprot (which are listed at <http://www.uniprot.org/faq/28#id_mapping_examples>), we can use `uniprot.py` to map the sequence identifiers to another type.

In this example, we have some RefSeq protein identifiers
(P_REFSEQ_AC) that we want to map to UniProt Accession
identifiers (ACC):

    seqids = "NP_000508.1  NP_001018081.3".split()

    pairs = uniprot.batch_uniprot_id_mapping_pairs(
      'P_REFSEQ_AC', 'ACC', seqids)

    pprint.pprint(pairs, indent=2)


## Getting protein sequence metadata

If we have our list of seqids in the Uniprot Accesion format, we can extract
the metadata of the proteins:

    uniprot_seqids = 'A0QSU3 D9QCH6 A0QL36'.split()
    uniprot_data = uniprot.batch_uniprot_metadata(
        uniprot_seqids, 'cache.basename')
    pprint.pprint(mapping, indent=2)

Now `batch_uniprot_metadata` uses a simple parser which only parses a small
number of fields in the metadata text, including the sequences for each
protein. For example, one common thing you can do is to save the sequences in
the metadata as a `.fasta` file. With the `uniprot_data` extracted from above,
you can write the output, using the convenience function:

    uniprot.write_fasta('output.fasta', uniprot_data, uniprot_seqids)

However, the raw uniprot metadata text is stored in files with the tag `cache.basename*.txt` and you can parse it yourself:

    for l in open('cache.basename0.txt'):
      print l

The information for the sequences of the isoforms for each protein are in the
metadta, but expressed in the form of VAR_SEQ entries. Here, we provide a
parser that reconstructs the isoform sequences directly from the metadata:
  
    text = open('metadata.cache0.txt').read()
    uniprot_data = uniprot.parse_isoforms(text)
    pprint.pprint(uniprot_data)

## Brute-force matching

Unfortunately, you probably have been given some files where you can't even recognize the type of sequence identifiers. You're won't to be able to get the  metadata, unless you can map your seqid to the Uniprot Accession id type:

Never fear, there's an brute-force method that you can use to figure out the id type of your seqids. Included is `seqidtype.py`, which should be installed as an executable script. On the command-line, just try:

    >> seqidtype.py YP_885981.1

`seqidtype.py` will run through all the seqid types listed in <http://www.uniprot.org/faq/28#id_mapping_examples> and try to get a uniprot accession mapping. Matches will be found only if the seqid type matches. Thus, after running through 50 odd individual matches, you will get a list of relevant seqid types. The output should look something like:

    ===> Analyzing YP_885981.1
    Fetching 1 (ACC->ACC) seqid mappings ...
    YP_885981.1 -> ACC -> None
    Fetching 1 (ID->ACC) seqid mappings ...
    YP_885981.1 -> ID -> None
    . 
    .
    .
    Fetching 1 (P_REFSEQ_AC->ACC) seqid mappings ...
    YP_885981.1:P_REFSEQ_AC -> A0QSU3
    .
    .
    .
    YP_885981.1 is compatible with seqid type: P_REFSEQ_AC

Since this requires a lot of internet calls, intermediate results are stored in the current directory under `seqidtype.json`. This may get big if you run the script a lot, and can be simply deleted. Once you have the seqid type in hand, you can now map your seqids to the Uniprot Accession id type:

    pairs = uniprot.batch_uniprot_id_mapping_pairs(
      'P_REFSEQ_AC', 'ACC', seqids)

## Filtering undecipherable seqid's

get_metadata_with_some_seqid_conversions(seqids, cache_fname=None)


As a bioinformatician, you are probably given FASTA files with sequence identifiers
from all sorts of different places. From these, you might extract an:

     seqids = open('seqids.txt', 'r').split()

However, most of the time we can cort the seqids into recognizable id types and then map each group appropriately. In `uniprot.py`, there is an example of a function that combines a mapping of several different types, before fetching thre uniprot metadata. 

## Chaining calls

Let's say you have a bunch of seqids of all different types. By chaining a bunch of calls to `uniprot.py`, you can construct a function that can map, store, and fetch metadata all in one go. An example is an included function that can take a bunch of seqids, identify and separate ENSEMBL, REFSEQ and UNIPROT ids, map to UNIPROT accession ids, and then fetch the metadata:

    metadata = uniprot.get_metadata_with_some_seqid_conversions(
         seqids, 'cache')

The heart of the function `get_metadata_with_some_seqid_conversions` uses pattern matching functions such as to `is_ensembl identify ENSEMBL ids, as can be seen in this fragment:

    # convert a few types into uniprot_ids
    id_types = [
      (is_sgd, 'ordered-locus-tag', 'ENSEMBLGENOME_PRO_ID'),
      (is_refseq, 'refseq', 'P_REFSEQ_AC'),
      (is_refseq, 'refseq', 'REFSEQ_NT_ID'),
      (is_ensembl, 'ensembl', 'ENSEMBL_ID'),
      (is_maybe_uniprot_id, 'uniprot-entity', 'ID')]
    for is_id_fn, name, uniprot_mapping_type in id_types:
      probe_id_type(entries, is_id_fn, name, uniprot_mapping_type, cache_fname+'.'+name)

The metadata is returned as a dictionary with the original seqids as keys.


(c) 2012, Bosco Ho


