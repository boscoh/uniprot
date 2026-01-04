

# UNIPROT.PY


`uniprot.py` provides a Python interface to the UniProt website <http://uniprot.org> to: 

1. Map between types of protein seqids (sequence identifiers)

2. Fetch metadata for proteins such as organism, sequence and GO annotations.


## Installation

`uniprot` is available as a PyPI library:

```bash
pip install uniprot
```

or

```bash
uv add uniprot
```

If you want to run from a local copy, download the package and sync dependencies:

```bash
uv sync
uv run python your_script.py
```

## Examples

Before you run any of the examples, import the module:

```python
import uniprot
```

It's super useful to import the pprint function to interrogate the data structures that the functions are returning:

```python
import pprint
```

A convenience function is provided to read seqids and sequences from a fasta file:

```python
seqids, fastas = uniprot.read_fasta('examples/example.fasta')
```


### Fetch seqid mappings

UniProt.org provides a seqid mapping service, but you must specify the seqid types, which are listed at <https://www.uniprot.org/help/id_mapping>. In this example, we have some RefSeq seqid's (RefSeq_Protein) that we want to map to UniProtKB identifiers:

```python
seqids = "NP_000508.1  NP_001018081.3".split()

pairs = uniprot.batch_uniprot_id_mapping_pairs(
  'RefSeq_Protein', 'UniProtKB', seqids)

pprint.pprint(pairs, indent=2)
```


### Getting protein sequence metadata

To get metadata for sequences, we need to have a list of seqids in the Uniprot Accesion or Uniprot ID format. To get the metadata:

```python
uniprot_seqids = 'A0QSU3 D9QCH6 A0QL36'.split()
uniprot_data = uniprot.batch_uniprot_metadata(
    uniprot_seqids, 'cache')
pprint.pprint(mapping, indent=2)
```

The function `batch_uniprot_metadata` contains a simple parser that extracts a small number of fields into a Python dictionary, with the Uniprot ID as the dictionary key. The results are obtained though batched queries to http://uniprot.org over several calls. An optional directory `cache` refers to a directory that stores cached results in case of interruption. You can carry further analysis on the `uniprot_data` dictionary. For example, you can write the sequences to a `.fasta` file using the convenience
function:

```python
uniprot.write_fasta('output.fasta', uniprot_data, uniprot_seqids)
```

If you would rather parse the metadata text yourself, you can refer to the raw text that was cached in the `cache/metadata.*.txt` files:

```python
for l in open('cache/metadata.0.txt'):
  print l
```

### Sorting seqids to find a good representative

Sometimes you will have a bunch of seqids that are related. For further analysis, you might just want pick the best one with the most useful uniprot information - for instance, the one that is the longest and that has also been reviewed (manually curated). 

A function `sort_seqids_by_uniprot` does just that. Let's say we have `uniprot_seqids` and `uniprot_data` from before. Then to find the most useful representative:

```python
sorted_seqids = uniprot.sort_seqids_by_uniprot(uniprot_seqids, uniprot_data)
best_seqid = sorted_seqids[0]
```

### Extracting isoform sequences

The Uniprot metadata contains information for the known isoforms of a protein, but this is expressed rather awkwardly as VAR_SEQ entries. Here is a function that reconstructs the isoform sequences from the raw metadata text:

```python
text = open('cache/metadata.0.txt').read()
isoforms_dict = uniprot.parse_isoforms(text)
pprint.pprint(isoforms_dict)
```


### Brute-force seqid-type matching

Unfortunately, you probably have been given some files where you can't recognize the seqid type. You are not going to be able to fetch the metadata unless you can map your seqid to the Uniprot Accession type.

Never fear! The `seqidtype_analyze()` function uses a brute-force approach to figure out the id type of a bunch of seqids. You can use it programmatically:

```python
uniprot.seqidtype_analyze('YP_885981.1', cache_fname='seqidtype.json')
```

Or run it as a command-line tool:

```bash
uv run seqidtype YP_885981.1
```

`seqidtype_analyze()` will attempt to map a seqid against all the seqid types listed in <https://www.uniprot.org/help/id_mapping>. After running through all ~100 seqid types, you will get a list of working seqid types, which should look something like:

```
Analyzing YP_885981.1
YP_885981.1:UniProtKB -> None
YP_885981.1:UniProtKB_AC-ID -> None
YP_885981.1:RefSeq_Protein -> A0QSU3
YP_885981.1 is compatible with: RefSeq_Protein
```

Since this requires lots of http requests, to avoid lost work, the intermediate results are cached in the current directory under `seqidtype.json`, which can be safely deleted. Once you have obtained the seqid type, you can map your seqids to the UniProtKB seqid type:

```python
pairs = uniprot.batch_uniprot_id_mapping_pairs(
  'RefSeq_Protein', 'UniProtKB', seqids)
```

## Chaining calls

Let's say you have a bunch of seqids of several different types. By chaining a bunch of calls to `uniprot.py`, you can construct a master function that fetches metadata for your seqids all in one go. Included is a function that can fetch metadata for ENSEMBL, REFSEQ and UNIPROT seqids:

```python
metadata = uniprot.get_metadata_with_some_seqid_conversions(
     seqids, 'cache')
```

The heart of the function `get_metadata_with_some_seqid_conversions` uses pattern matching functions, such as `is_ensembl` to identify ENSEMBL ids, as can be seen in this fragment:

```python
id_types = [
  (is_sgd, 'locustag', 'SGD'),
  (is_refseq, 'refseqp', 'RefSeq_Protein'),
  (is_refseq, 'refseqnt', 'RefSeq_Nucleotide'),
  (is_ensembl, 'ensembl', 'Ensembl'),
  (is_maybe_uniprot_id, 'uniprotid', 'UniProtKB_AC-ID')]
for is_id_fn, name, uniprot_mapping_type in id_types:
  probe_id_type(entries, is_id_fn, name, uniprot_mapping_type, cache_fname+'.'+name)
```

The metadata is then returned as a dictionary with the original seqids as keys. You can follow the logic in this function to construct functions of your own design.

## Project Structure

```
uniprot/
├── uniprot.py          # Main module
├── pyproject.toml      # Project configuration
├── readme.md
├── tests/              # Test suite
│   ├── test_uniprot.py       # Unit tests
│   ├── test_integration.py   # API integration tests
│   └── data/                 # Test fixtures
│       └── isoform/          # Isoform test data
└── examples/           # Usage examples
    ├── example.py
    └── example.fasta
```

- **uniprot.py** - Main module with functions for ID mapping, metadata fetching, and parsing
- **tests/test_uniprot.py** - Unit tests for parsing, FASTA I/O, and ID detection
- **tests/test_integration.py** - Integration tests with real UniProt API (requires network)
- **examples/** - Example usage demonstrating the module's functionality

## Testing

Run the unit test suite with uv:

```bash
uv run python -m unittest tests.test_uniprot -v
```

For integration tests (requires internet):

```bash
uv run python -m unittest tests.test_integration -v
```

Run all tests:

```bash
uv run python -m unittest discover tests -v
```

## Changelog

### 1.4.1 (January 4, 2026)
- Reorganized project structure: tests in `tests/`, examples in `examples/`
- Fixed bug in `get_metadata_with_some_seqid_conversions()` where empty seqids caused HTTP 400 errors
- Fixed VAR_SEQ parsing to handle both old and new UniProt formats
- Fixed isoform parsing for entries without `(in isoform X)` annotations
- Fixed parsing of VAR_SEQ entries without `->` transitions

### 1.4 (January 4, 2026)
- Uses `pyproject.toml` for project configuration
- Dependency management with [uv](https://docs.astral.sh/uv/)
- Moved seqidtype functionality into uniprot module as `seqidtype_analyze()` and `seqidtype_cli()`
- **Test Suite**: Comprehensive test coverage
  - `test_uniprot.py` - Unit tests for parsing, FASTA I/O, and ID detection
  - `test_integration.py` - Integration tests with real UniProt API endpoints
- **API Update**: Updated to use new UniProt REST API field names (July 2021+)
  - Old field names (e.g., `P_REFSEQ_AC`, `ENSEMBL_ID`, `ACC`, `ID`) are deprecated
  - New field names (e.g., `RefSeq_Protein`, `Ensembl`, `UniProtKB_AC-ID`, `UniProtKB`) are now used
  - Added dynamic validation: code now fetches mapping rules from the API to validate field combinations before requests
  - Maps to `UniProtKB` destination instead of `UniProtKB_AC-ID` (which is source-only)

### 1.3.2 (November 22, 2021)
- Update url for mapping (knightjdr)

### 1.3.1 (May 17, 2017)
- Added support for TrEMBL format uniprot accessions (Peter Oxley)

### 1.3 (February 1, 2016)
- Migrated to Python 3
- Added isoform parsing support
- Improved seqid type detection

### 1.2 (March 26, 2015)
- changed the cache parameter of `batch_uniprot_id_mapping_pairs` and `batch_uniprot_metadata`  to a directory `cache_dir`
- the batch functions now saves the seqids parameters and will do a clean search  if the cached seqids do not match
- abstracted all screen output to the `logging` function that can be overwritten

### 1.1 (February 18, 2015)
- sort_seqids_by_uniprot
- change limits to 400 (due to some error messages from uniprot)

### 1.0.2 (July 23, 2014)
- add a default cache_fname parameter to get_uniprot_id_mapping_pairs 

### 1.0.1 (February 26, 2014)
- bug parsing isoform metadata when dangling isoforms at the end of line
- get_metadata_with_some_seqid_conversions can now actually handle None for cache_basename

### 1.0.0 (February 7, 2014)
- Initial release

(c) 2013, Bosco Ho


