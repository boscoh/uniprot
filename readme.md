

# UNIPROT.PY


`uniprot.py` provides a Python interface to the UniProt website <http://uniprot.org> to: 

1. Map between types of protein seqids (sequence identifiers)

2. Fetch metadata for proteins such as organism, sequence and GO annotations.


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
        uniprot_seqids, 'cache')
    pprint.pprint(mapping, indent=2)

The function `batch_uniprot_metadata` contains a simple parser that extracts a small number of fields into a Python dictionary, with the Uniprot ID as the dictionary key. The results are obtained though batched queries to http://uniprot.org over several calls. An optional directory `cache` refers to a directory that stores cached results in case of interruption. You can carry further analysis on the `uniprot_data` dictionary. For example, you can write the sequences to a `.fasta` file using the convenience
function:

    uniprot.write_fasta('output.fasta', uniprot_data, uniprot_seqids)

If you would rather parse the metadata text yourself, you can refer to the raw text that was cached in the `cache/metadata.*.txt` files:

    for l in open('cache/metadata.0.txt'):
      print l

### Sorting seqids to find a good representative

Sometimes you will have a bunch of seqids that are related. For further analysis, you might just want pick the best one with the most useful uniprot information - for instance, the one that is the longest and that has also been reviewed (manually curated). 

A function `sort_seqids_by_uniprot` does just that. Let's say we have `uniprot_seqids` and `uniprot_data` from before. Then to find the most useful representative:
 
    sorted_seqids = uniprot.sort_seqids_by_uniprot(uniprot_seqids, uniprot_data)
    best_seqid = sorted_seqids[0]

### Extracting isoform sequences

The Uniprot metadata contains information for the known isoforms of a protein, but this is expressed rather awkwardly as VAR_SEQ entries. Here is a function that reconstructs the isoform sequences from the raw metadata text:
  
    text = open('cache/metadata.0.txt').read()
    isoforms_dict = uniprot.parse_isoforms(text)
    pprint.pprint(isoforms_dict)


### Brute-force seqid-type matching

Unfortunately, you probably have been given some files where you can't recognize the seqid type. You are not going to be able to fetch the metadata unless you can map your seqid to the Uniprot Accession type.

Never fear!  Included is `seqidtype`, an executable script that uses a  brute-force approach to figure out the id type of a bunch of seqids. On the command-line:

    >> seqidtype YP_885981.1

`seqidtype` will attempt to map a seqid against all the seqid types listed in <http://www.uniprot.org/faq/28#id_mapping_examples>. After running through all 50 or so seqid types, you will get a list of working seqid types, which should look something like:

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

The heart of the function `get_metadata_with_some_seqid_conversions` uses pattern matching functions, such as `is_ensembl` to identify ENSEMBL ids, as can be seen in this fragment:

    # convert a few types into uniprot_ids
    id_types = [
      (is_sgd, 'locustag', 'ENSEMBLGENOME_PRO_ID'),
      (is_refseq, 'refseqp', 'P_REFSEQ_AC'),
      (is_refseq, 'refseqnt', 'REFSEQ_NT_ID'),
      (is_ensembl, 'ensembl', 'ENSEMBL_ID'),
      (is_maybe_uniprot_id, 'uniprotid', 'ID')]
    for is_id_fn, name, uniprot_mapping_type in id_types:
      probe_id_type(entries, is_id_fn, name, uniprot_mapping_type, cache_fname+'.'+name)

The metadata is then returned as a dictionary with the original seqids as keys. You can follow the logic in this function to construct functions of your own design.

## Changelog

### 1.3
- Python 3 compatibility (thanks Nader Moshed)

### 1.2
- changed the cache parameter of `batch_uniprot_id_mapping_pairs` and `batch_uniprot_metadata`  to a directory `cache_dir`
- the batch functions now saves the seqids parameters and will do a clean search  if the cached seqids do not match
- abstracted all screen output to the `logging` function that can be overwritten

### 1.1
- sort_seqids_by_uniprot
- change limits to 400 (due to some error messages from uniprot)

### 1.0.2
- add a default cache_fname parameter to get_uniprot_id_mapping_pairs 

### 1.0.1 
- bug parsing isoform metadata when dangling isoforms at the end of line
- get_metadata_with_some_seqid_conversions can now actually handle None for cache_basename

(c) 2013, Bosco Ho


