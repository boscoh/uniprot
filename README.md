

UNIPROT.PY
==========

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



# Examples

Before you run any of the examples, import the module:

    import uniprot

I also like importing pprint to print out dictionaries,
in order to see what I am getting back:

    import pprint

Reading seqids from a fasta file:

    seqids, fastas = uniprot.read_fasta('example.fasta')


## Fetching identifier mappings

If we have a bunch of seqids that cleanly maps to a known 
type of id recognized by uniprot, we can use batch mode to
map identifiers in groups of 100 at a time. It's very important
that the type of the identifiers are specified. These can be
looked up at <http://www.uniprot.org/faq/28#mapping-faq-table>. 

In this example, we have some RefSeq protein identifiers
(P_REFSEQ_AC) that we want to map to UniProt Accession
identifiers (ACC):

    seqids = " NP_000508.1  NP_001018081.3 ".split()

    pairs = uniprot.batch_uniprot_id_mapping_pairs(
      'P_REFSEQ_AC', 'ACC', seqids)

    pprint.pprint(pairs, indent=2)

Most of the routines require a caching filename as parameter.
This is because such queries are very temperamental - it
depends on your network latency, as well as the good graces
of uniprot.org itself. As such, the functions caches
as much on disk as possible to avoid lost work.


## Brute-force matching

If we have a bunch of seqids that can't be recognized by
the batch identifier-mapping service of UniProt (happens way more
often than you'd think), you can use the sequential service
which is much more forgiving. It does not need a prespecification
of the identifier type. However, the tradeoff is that
only one identifier can be looked up at a time, so it is
excruciatingly slow:

    seqids = """
    EFG_MYCA1 YP_885981.1 CAS1_BOVIN CpC231_1796
    """.split()

    mapping = uniprot.sequentially_convert_to_uniprot_id(
        seqids, 'cache.dict')

    pprint.pprint(mapping, indent=2)

## Getting protein sequence metadata

If we have our list of uniprot seqids, we can then extract
the metadata, which includes, for instance the sequence of the
protein:

    uniprot_seqids = 'A0QSU3 D9QCH6 A0QL36'.split()
    uniprot_data = uniprot.batch_uniprot_metadata(
        uniprot_seqids, 'cache.txt')
    pprint.pprint(mapping, indent=2)

Now `batch_uniprot_metadata` uses a simple uniprot text parser
which only parses a small number of fields. 

Of course, you'd probably want to weigh in with your UniProt parser.
Simply read in the cache.txt text file:

    for l in open('cache.txt'):
      print l

For instance, you might want to save the fasta sequences:

    uniprot.write_fasta('output.fasta', uniprot_data, uniprot_seqids)

## Chaining calls

By chaining a bunch of calls, you can construct functions
that directly transform ID's of your choice:

    import os

    def map_to_refseq(seqids):
      uniprot_mapping = uniprot.sequentially_convert_to_uniprot_id(
          seqids, 'func.cache.dict')
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
      os.remove('func.cache.dict')
      return mapping

    seqids = " EFG_MYCA1 YP_885981.1 CpC231_1796 ".split()
    mapping = map_to_refseq(seqids)
    pprint.pprint(mapping, indent=2)



(c) 2012, Bosco Ho


