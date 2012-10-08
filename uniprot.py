import os
import json
import textwrap
import time
import StringIO
import requests

"""
This is my uniprot python library. It provides a bunch of 
routines to access the uniprot.org website's RESTful API.

Principally, you want to do two things:
1) Map different types of sequence identifiers
2) Fetch uniprot metadata for certain sequence id's.

Given a set of seqid's, this library allows a programmic
approach to transforming id's and if the uniprot_id is
obtained, it allows the extraction of the relevant metadata,
including the protein sequence.
"""


# Doing internet lookups is fraught with error, 
# so caching to disk is recommended. Here's a
# couple of helper functions to cache Python dictionaries.


def read_dict(cache_json):
  return json.load(open(cache_json, 'rU'))


def write_dict(a_dict, cache_json, indent=2):
  json.dump(a_dict, fp=open(cache_json, 'w'), indent=indent)


def cache_dict_fn(dict_fn, cache_json, is_overwrite=False):
  """
  This is a convenience wrapper for lazy evaluation
  of a function that returns a dictionary. If the
  cache_json exists, the file is read and returned,
  otherwise the fn() is run.
  """
  if not is_overwrite and os.path.isfile(cache_json):
    return read_dict(cache_json)
  else:
    results = dict_fn()
    write_dict(results, cache_json)
    return results



def batch_uniprot_id_mapping_pairs(
    from_type, to_type, seqids, n_batch=100):
  """
  Converts identifiers using above function, but launches url
  requests in safe batches of 100 seqids at a time.
  
  from_type and to_type can be obtained from:
    http://www.uniprot.org/faq/28#mapping-faq-table

  Returns a list of matched pairs of identifiers.
  """
  r = requests.post(
      'http://www.uniprot.org/mapping/', 
      params={
        'from': from_type,
        'to': to_type,
        'format': 'tab',
        'query': ' '.join(seqids)})
  lines = [l for l in r.text.splitlines() if 'from' not in l.lower()]
  return [l.split('\t')[:2] for l in lines]



def sequentially_convert_to_uniprot_id(seqids, cache_json):
  """
  This is the key function to convert an arbitary seqid to a
  uniprot id. It only works for one seqid per url request, so
  it is slow. But if you cannot transform your seqids so that
  batch mapping can work, it is the only option. 

  Since this is slow, it caches the results at every pass.

  Returns a dictionary of the input seqids as keys. 
  """
  if not os.path.isfile(cache_json):
    mapping = {}
  else:
    mapping = read_dict(cache_json)
  for l in seqids:
    seqid = l.split()[0]
    if seqid not in mapping:
      result = batch_uniprot_id_mapping_pairs(
            'ACC+ID', 'ACC', ['?' + seqid])
      if result:
        mapping[seqid] = result[0][0]
        write_dict(mapping, cache_json)
        print seqid, "->", mapping[seqid]
      else:
        print seqid, '-> [null]'
  return mapping



def parse_uniprot_txt_file(cache_txt):
  """
  Parses a text file of metadata retrieved from uniprot.org.

  Only a few fields have been parsed, but this provides a
  template for the other fields.

  A single description is generated from joining alternative
  descriptions.

  Returns a dictionary with the uniprot id as keys.
  """

  tag = None
  uniprot_id = None
  uniprot_data = {}
  for l in open(cache_txt):
    test_tag = l[:5].strip()
    if test_tag and test_tag != tag:
      tag = test_tag
    line = l[5:].strip()
    words = line.split()
    if tag == "ID":
      uniprot_id = words[0]
      is_reviewed = words[1].startswith('Reviewed')
      length = int(words[2])
      uniprot_data[uniprot_id] = {
        'id': uniprot_id,
        'is_reviewed': is_reviewed,
        'length': length,
        'sequence': '',
        'accs': [],
      }
      entry = uniprot_data[uniprot_id]
    if tag == "SQ":
      if words[0] != "SEQUENCE":
        entry['sequence'] += ''.join(words)
    if tag == "AC":
      accs = [w.replace(";", "") for w in words]
      entry['accs'].extend(accs)
    if tag == "DR":
      if 'RefSeq' in words[0]:
        if 'refseq' not in entry:
          entry['refseq'] = []
        ids = [w[:-1] for w in words[1:]]
        entry['refseq'].extend(ids)
      if 'KEGG' in words[0]:
        if 'kegg' not in entry:
          entry['kegg'] = []
        ids = [w[:-1] for w in words[1:]]
        ids = filter(lambda w: len(w) > 1, ids)
        entry['kegg'].extend(ids)
    if tag == "GN":
      if 'gene' not in entry and len(words) > 0:
        pieces = words[0].split("=")
        if len(pieces) > 1 and 'name' in pieces[0].lower():
          entry['gene'] = pieces[1].replace(';', '').replace(',', '')
    if tag == "OS":
      entry['organism'] = line
    if tag == "DE":
      if 'descriptions' not in entry:
        entry['descriptions'] = []
      entry['descriptions'].append(line)

  for entry in uniprot_data.values():
    descriptions = entry['descriptions']
    for i in reversed(range(len(descriptions))):
      description = descriptions[i]
      if 'Short' in description or 'Full' in description:
        j = description.find('=')
        descriptions[i] = description[j+1:].replace(';', '')
      else:
        del descriptions[i]
    entry['description'] = '; '.join(descriptions)

  return uniprot_data



def batch_uniprot_metadata(seqids, cache_txt):
  """
  Returns a dictonary of the uniprot metadata (as parsed 
  by parse_uniprot_txt_file) of the given seqids. The seqids
  must be valid uniprot identifiers.

  Returns a dictionary of uniprot data with the input
  seqids as keys.
  """

  if os.path.isfile(cache_txt):
    print "Loading uniprot data from previous lookup", cache_txt
  else:
    print "Looking up %d seqids at http://www.uniprot.org" % len(seqids)
    r = requests.post(
        'http://www.uniprot.org/batch/', 
        files={'file':StringIO.StringIO(' '.join(seqids))}, 
        params={'format':'txt'})
    while 'Retry-After' in r.headers:
      t = int(r.headers['Retry-After'])
      print 'Waiting %d' % t
      time.sleep(t)
      r = requests.get(r.url)
    open(cache_txt, 'w').write(r.text)

  print "Reading from", cache_txt
  uniprot = parse_uniprot_txt_file(cache_txt)

  # resort the dictionary wrt to input seqids as keys
  results = {}
  for uniprot_id in uniprot.keys():
    for seqid in uniprot[uniprot_id]['accs']:
      results[seqid] = uniprot[uniprot_id]
  for seqid in seqids:
    if seqid not in results:
      primary_seqid = seqid[:6]
      if primary_seqid in results:
        results[seqid] = results[primary_seqid]
  return results



def get_filtered_uniprot_metadata(seqids, cache_txt):
  """
  Returns a dictionary of uniprot data, but filters
  seqids for uniprot identifiers by using a mapping call
  to uniprot first.
  """

  stripped_seqids = [s[:6] for s in seqids]
  pairs = batch_uniprot_id_mapping_pairs(
      'ACC', 'ACC', stripped_seqids)
  uniprot_seqids = []
  for seqid1, seqid2 in pairs:
    if seqid1 in stripped_seqids and seqid1 not in uniprot_seqids:
      uniprot_seqids.append(seqid1)
  uniprot_dict = batch_uniprot_metadata(uniprot_seqids, cache_txt)
  for seqid in seqids:
    if seqid not in uniprot_seqids and seqid[:6] in uniprot_seqids:
      uniprot_dict[seqid] = uniprot_dict[seqid[:6]]
  return uniprot_dict



def sort_seqids_by_uniprot(names, uniprot_dict):
  """
  Weighs protein names to whether they are found
  by uniprot, and are reviewed
  """
  def seqid_val(seqid):
    val = 3
    if seqid in uniprot_dict:
      val -= 1
      if uniprot_dict[seqid]['is_reviewed']:
        val -= 1
      if len(seqid) <= 6:
        val -= 1
    return val
  names.sort(cmp=lambda a, b: seqid_val(a) - seqid_val(b))
  return names    


def parse_fasta_header(header, seqid_fn=None):
  """
  Parses a FASTA format header (with our without the initial '>').
  
  If NCBI SeqID format (gi|gi-number|gb|accession etc, is detected
  the first id in the list is used as the canonical id (see see
  http://www.ncbi.nlm.nih.gov/books/NBK21097/#A631 ).

  Extra processing to parse the seqid can be provided
  by giving a seqid_fn function.

  Returns a tuple of sequence id and sequence name/description.
  """
  # check to see if we have an NCBI-style header
  if header[0] == '>':
    header = header[1:]
  if header.find("|") != -1:
    tokens = header.split('|')
    # "gi|ginumber|gb|accession bla bla" becomes "gi|ginumber"
    seqid = "%s|%s" % (tokens[0], tokens[1].split()[0])
    name = seqid + ' ' + tokens[-1:][0].strip()
  else:
    # otherwise just split on spaces & hope for the best
    tokens = header.split()
    seqid = tokens[0]
    name = header[0:-1].strip()
  
  if seqid_fn is not None:
    seqid = seqid_fn(seqid)

  return seqid, name



def read_selected_fasta(seqids, fasta_db, seqid_fn=None):
  """
  Extracts protein sequences from a fasta database, given
  a list of desired seqids. A seqid_fn can be given that
  parses both the input seqids and the seqid in the fasta 
  file to faciliate matching.

  Returns a dictionary for each seqid where the sequence
  is found in the 'sequence' field.
  """
  live_name = None
  proteins = {}
  if seqid_fn is not None:
    seqids = map(seqid_fn, seqids)
  for i, line in enumerate(open(fasta_db)):
    if line.startswith(">"):
      fasta_seqid, description = \
          parse_fasta_header(line, seqid_fn)
      live_name = None
      for seqid in seqids:
        if fasta_seqid == seqid:
          live_name = fasta_seqid
          proteins[live_name] = {
            'sequence': "",
            'description': description,
          }
          break
    elif live_name:
      proteins[live_name]['sequence'] += line.strip()
  for seqid in proteins:
    sequence = proteins[seqid]['sequence']
    if sequence:
      proteins[seqid]['length'] = len(sequence)
  return proteins



def read_fasta(fasta_db, seqid_fn=None):
  """
  Parses a fasta database file. Warning: will be slow for very large
  databases. In that case, it would be better to use 
  read_selected_fasta() instead.

  Returns a lsit of seqids encountered, and a dictionary
  of sequences for each seqid.
  """
  seqids = []
  seqid = None
  proteins = {}
  for line in open(fasta_db):
    if line.startswith(">"):
      seqid, description = parse_fasta_header(line, seqid_fn)
      seqids.append(seqid)
      proteins[seqid] = {
        'sequence':"",
        'description':description,
      }
      continue
    if seqid is not None:
      words = line.split()
      if words:
        proteins[seqid]['sequence'] += words[0]
  return seqids, proteins
  


def write_fasta(
    fasta_fname, proteins, seqids, width=50):
  """
  Creates a fasta file of the sequences of a subset of the proteins.
  """
  f = open(fasta_fname, "w")
  for seqid in seqids:
    seq_wrap = textwrap.fill(proteins[seqid]['sequence'], width)
    f.write(">%s %s\n%s\n" % (seqid, proteins[seqid]['description'], seq_wrap))
  f.close()



