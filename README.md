

# UNIPROT.PY


`uniprot.py` provides a Python interface to the UniProt website <http://uniprot.org> to: 

1. Map between types of protein seqids (sequence identifiers)

2. Fetch metadata for proteins such as organism, the sequence and GO annotations.


## Installation

To install:

    >> pip install uniprot

If you're installing directly from `setup.py`, be aware that `uniprot.py` has one dependency: `requests.py`:

    >> pip install requests

## Examples

Before you run any of the examples, import the module:

    import uniprot

It's super useful to import the pprint function to interrogate the data structures that the functions are returning:

    import pprint

A convenience function is provided to read seqids and sequences from a fasta file:

    seqids, fastas = uniprot.read_fasta('example.fasta')


### Fetch seqid mappings

Uniprot.org provides a seqid mapping service, but you must specify the seqid types, which are listed at <http://www.uniprot.org/faq/28#id_mapping_examples>.  In this example, we have some RefSeq seqid's (P_REFSEQ_AC) that we want to map to UniProt Accession seqid's (ACC):

    seqids = "NP_000508.1  NP_001018081.3".split()

    pairs = uniprot.batch_uniprot_id_mapping_pairs(
      'P_REFSEQ_AC', 'ACC', seqids)

    pprint.pprint(pairs, indent=2)


### Getting protein sequence metadata

To get metadata for sequences, we need to have a list of seqids in the Uniprot Accesion or Uniprot ID format. To get the metadata:

    uniprot_seqids = 'A0QSU3 D9QCH6 A0QL36'.split()
    uniprot_data = uniprot.batch_uniprot_metadata(
        uniprot_seqids, 'cache.basename')
    pprint.pprint(mapping, indent=2)

The function `batch_uniprot_metadata` contains a simple parser that extracts a small number of fields into a Python dictionary, with the Uniprot ID as the dictionary key. You can carry further analysis on this dictionary. For example, you can write the sequences to a `.fasta` file using the convenience
function:

    uniprot.write_fasta('output.fasta', uniprot_data, uniprot_seqids)

The metadata contains information for the known isoforms of a protein, but this is expressed rather awkwardly as VAR_SEQ entries. Here is a function that reconstructs the isoform sequences from the metadata:
  
    text = open('metadata.cache0.txt').read()
    isoforms_dict = uniprot.parse_isoforms(text)
    pprint.pprint(isoforms_dict)

If you would rather parse the metadata text yourself, the raw text is cached
to `cache.basename*.txt` files:

    for l in open('cache.basename0.txt'):
      print l


### Brute-force seqid-type matching

Unfortunately, you probably have been given some files where you can't recognize the seqid type. You are not going to be able to fetch the metadata unless you can map your seqid to the Uniprot Accession type.

Never fear!  Included is `seqidtype.py`, an executable script that uses a  brute-force method approach to figure out the id type of a bunch of seqids. On the command-line:

    >> seqidtype.py YP_885981.1

`seqidtype.py` will attempt to map a seqid against all the seqid types listed in <http://www.uniprot.org/faq/28#id_mapping_examples>. After running through all 50 or so seqid types, you will get a list of working seqid types, which should look something like:

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

Since this requires lots of http requests, to avoid lost work, the intermediate results are cached in the current directory under `seqidtype.json`, which can be safely deleted. Once you have obtained the seqid type, you can map your seqids to the Uniprot Accession seqid type:

    pairs = uniprot.batch_uniprot_id_mapping_pairs(
      'P_REFSEQ_AC', 'ACC', seqids)

## Chaining calls

Let's say you have a bunch of seqids of several different types. By chaining a bunch of calls to `uniprot.py`, you can construct a master function that fetches metadata for your seqids all in one go. Included is a function that can fetch metadata for ENSEMBL, REFSEQ and UNIPROT seqids:

    metadata = uniprot.get_metadata_with_some_seqid_conversions(
         seqids, 'cache')

The heart of the function `get_metadata_with_some_seqid_conversions` uses pattern matching functions, such as to `is_ensembl` to identify ENSEMBL ids, as can be seen in this fragment:

    # convert a few types into uniprot_ids
    id_types = [
      (is_sgd, 'ordered-locus-tag', 'ENSEMBLGENOME_PRO_ID'),
      (is_refseq, 'refseq', 'P_REFSEQ_AC'),
      (is_refseq, 'refseq', 'REFSEQ_NT_ID'),
      (is_ensembl, 'ensembl', 'ENSEMBL_ID'),
      (is_maybe_uniprot_id, 'uniprot-entity', 'ID')]
    for is_id_fn, name, uniprot_mapping_type in id_types:
      probe_id_type(entries, is_id_fn, name, uniprot_mapping_type, cache_fname+'.'+name)

The metadata is then returned as a dictionary with the original seqids as keys. You can follow the logic in this function to construct functions of your own design.

(c) 2013, Bosco Ho


