

# UNIPROT.PY


This is my python library for dealing with sequence identifier
mappings, which leverages the UniProt website http://uniprot.org,
which is arguably the best protein sequence resource out there.

This library provides a interface to access the 
RESTful API of the uniprot.org website.

Principally, you want to do two things:

1. Map between different types of sequence identifiers

2. Fetch uniprot metadata for protein sequence identifiers, 
including organism, the protein sequence itself, GO annotations.

This library provides a simple UniProt parser to analyze the
sequence metadata. It is easy to switch in your own parser.

There are various ways to access the mapping service, 
where each method has its pros and cons (speed vs. robustness).


## Installation

The library `uniprot.py` only has one dependency, the wonderful requests.py library which wraps the truly dire httplib in standard Python library. To install:

    >> pip install uniprot

And if just to make sure:
    
    >> pip install requests 


## Examples

Before you run any of the examples, import the module:

    import uniprot

I also like importing pprint to print out dictionaries,
in order to see what I am getting back:

    import pprint

Reading seqids from a fasta file:

    seqids, fastas = uniprot.read_fasta('example.fasta')


## Fetching identifier mappings

If we have a bunch of sequence identifiers that cleanly maps to a known 
identifier type recognized by UniProt, we can use batch mode to
map identifiers in groups of 100 at a time. It's very important
that the type of the identifiers are specified. These can be
looked up at <http://www.uniprot.org/faq/28#id_mapping_examples>. 

In this example, we have some RefSeq protein identifiers
(P_REFSEQ_AC) that we want to map to UniProt Accession
identifiers (ACC):

    seqids = " NP_000508.1  NP_001018081.3 ".split()

    pairs = uniprot.batch_uniprot_id_mapping_pairs(
      'P_REFSEQ_AC', 'ACC', seqids)

    pprint.pprint(pairs, indent=2)

## Brute-force matching

Let's say you a bunch of sequence identifiers that you can't recognize (happens quite a lot, no?), then you're not going to be able to get the correct uniprot accession id's to fetch uniprot metadata :(.

Never fear, there's an ugly brute-force method that you can use to figure out what uniprot sequence-id type that your seqids have. Included in the package is a python script called `seqidtype.py`, which should be installed as an executable script, on the command-line, just try:

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

Since this is a laborious internet process, intermediate results are stored in the current directory under `seqidtype.json`. This may get big if you run the script a lot, and can be simply deleted. Once you have the seqid type in hand, you can now rerun the batch seqid mapping method.


## Getting protein sequence metadata

If we have our list of seqids in uniprot acc format, we can extract
the metadata of the proteins, including the sequence:protein:

    uniprot_seqids = 'A0QSU3 D9QCH6 A0QL36'.split()
    uniprot_data = uniprot.batch_uniprot_metadata(
        uniprot_seqids, 'cache.basename')
    pprint.pprint(mapping, indent=2)

Now `batch_uniprot_metadata` uses a simple uniprot text parser
which only parses a small number of fields. The actual uniprot metadata
is stored in files that begin with 'cache.basename*.txt', depending
on how many different batches are required.

Of course, you'd probably want to weigh in with your UniProt parser. If 
there is only one batch done, then only one file needs to be read, 
`cache.basename0.txt'. Simply read in this file:

    for l in open('cache.basename0.txt'):
      print l

One common thing you can do is to save the sequences in the metadata
as a .fasta file. With the `uniprot_data` extracted from above, you
can write the output as:

    uniprot.write_fasta('output.fasta', uniprot_data, uniprot_seqids)


## Filtering undecipherable seqid's

get_metadata_with_some_seqid_conversions(seqids, cache_fname=None)


As a bioinformatician, you are probably given FASTA files with sequence identifiers
from all sorts of different places. From these, you might extract an:

     seqids = open('seqids.txt', 'r').split()

Looking up the protein metadata using uniprot is difficult because non-uniprot
ids could produce different kinds of behavior. The best way is simply to
call a seemingly redundant call of ID mapping from uniprot to uniprot:

    pairs = uniprot.batch_uniprot_id_mapping_pairs('ACC+ID', 'ACC', seqids)

This will pick out the seqids that uniprot.org can recognize. This is a 
common enough operation to warrant its own function 
  
    metadata = uniprot.get_filtered_uniprot_metadata(seqids)

where a preemptive call to the uniprot mapping is used to filter out non-uniprot
identifiers.

Finally, an example of a more comprehensive metadata function is provided that
does some limited parsing for some more common ID types. Currrently, pattern
matching is used to identify ENSEMBL, REFSEQ and YEAST ORF seqids, as well,
the function can handle the typical xx|xxxxxxx|xxx format found in many databases.
Simply run the seqids through 

    metadata = uniprot.get_metadata_with_some_seqid_conversions(
         seqids, 'cache.txt')


## Chaining calls

By chaining a bunch of calls, you can construct functions
that directly transform ID's of your choice:

    import os

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

    seqids = " EFG_MYCA1 YP_885981.1 CpC231_1796 ".split()
    mapping = map_to_refseq(seqids)
    pprint.pprint(mapping, indent=2)



(c) 2012, Bosco Ho


