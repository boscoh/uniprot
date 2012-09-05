

UNIPROT.PY
==========

This is my UNIPROT python library for dealing with 
protein sequence identifiers. It leverages uniprot.org,
arguably the best protein sequence resource out there.

This library provides a interface to access the 
RESTful API of the uniprot.org website.

Principally, you want to do two things:

1. Map between different types of sequence identifiers

2. Fetch uniprot metadata for protein sequence identifiers, 
including organism, the protein sequence itself, GO annotations.

This library provides a simple uniprot parser, which you
can easily add to, or modify.

There are various ways to access the mapping service, 
which depends on the limitations of the service.



# Examples

Before you run any of the examples, import the module:

    import uniprot

And I also like importing pprint to print out dictionaries,
in order to see what I am getting back:

    import pprint

Reading seqids from a fasta file:

    seqids, fastas = uniprot.read_fasta('c_gl.fasta')

## Fetching identifier mappings

If we have a bunch of seqids that cleanly maps to a known 
type of id recognized by uniprot, we can use batch mode to
map id's in groups of 100 at a time:

    seqids = """
    
    """.split()
    mapping = uniprot.batch_uniprot_id_mapping_pairs(
      'ENSEMBL', 'ACC', seqids, 'cache.mapping.dict')

    pprint.pprint(mapping, indent=2)

Most of the routines require a caching filename as parameter.
This is because such queries are very temperamental - it
depends on your network latency, as well as the good graces
of uniprot.org itself. As such, the routines here caches
as much on disk as possible to avoid lost work.


## Brute-force matching

If we have a bunch of seqids that can't be recognized by
the batch mapping ID service of uniprot (happens way more
often than you'd think), you can use the sequential service
which is much more forgiving, but can only spit out uniprot
sequence identifiers:

    seqids = """
    
    """.split()
    mapping = uniprot.sequentially_convert_to_uniprot_id(
        seqids, fname_dict)
    pprint.pprint(mapping, indent=2)

## Getting protein sequence metadata

If we have our list of uniprot seqids, we can then extract
the metadata, which includes, for instance the sequence of the
protein:

    uniprot_data = uniprot.batch_uniprot_metadata(
        uniprot_seqids, 'c_gl.uniprot.txt')

Now `batch_uniprot_metadata` uses a simple uniprot text parser
which only parses a small number of fields. If you want to
weigh in with your uniprot parser, then you can roll your
own parser:

    uniprot_data = uniprot.batch_uniprot_metadata_to_txt_file(
        uniprot_seqids, 'cache.txt')
    for l in open('cache.txt'):
      print l

## Chaining calls

By chaining a bunch of calls, you can construct functions
that directly transform ID's of your choice:

    def get_mapping(seqids):
      mapping1 = uniprot.sequentially_convert_to_uniprot_id(
          seqids, 'cache1.dict')
      uniprot_ids = seqids.values()
      mapping2 = uniprot.batch_uniprot_id_mapping_pairs(
        'ACC', 'KEGG', uniprot_ids, 'cache2.dict')
      return mapping


(c) 2012, Bosco Ho


