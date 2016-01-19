# Copyright (c) 2013, Bosco Ho. All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import re
import sys
import os
import textwrap
import time
import json
try:
  from StringIO import StringIO
except ImportError:
  from io import StringIO
import shutil
from copy import deepcopy

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


def logging(s):
  sys.stdout.write(s)


def is_html(text):
  if re.search('<html', text):
    return True
  return False


def get_uniprot_id_mapping_pairs(
    from_type, to_type, seqids, cache_fname=None):
  """
  Returns a list of matched pairs of identifiers.

  from_type and to_type can be obtained from:
    http://www.uniprot.org/faq/28#mapping-faq-table
  """
  if cache_fname and os.path.isfile(cache_fname):
    logging("Loading cached (%s->%s) mappings from %s\n" % (from_type.upper(), to_type.upper(), cache_fname))
    text = open(cache_fname).read()
  else:
    logging("Fetching %s (%s->%s) mappings from http://uniprot.org...\n" % (len(seqids), from_type.upper(), to_type.upper()))
    r = requests.post(
        'http://www.uniprot.org/mapping/',
         files={'file':StringIO(' '.join(seqids))},
         params={
          'from': from_type.upper(),
          'to': to_type.upper(),
          'format': 'tab',
          'query': ''})
    text = r.text
    if cache_fname:
      with open(cache_fname, 'w') as f:
        f.write(text)
  if is_html(text):
    # failed call results in a HTML error reporting page
    logging("Error in fetching metadata\n")
    return []
  lines = [l for l in text.splitlines() if 'from' not in l.lower()]
  return [l.split('\t')[:2] for l in lines]


def batch_uniprot_id_mapping_pairs(
    from_type, to_type, seqids, batch_size=400, cache_dir=None):
  """
  Returns a list of matched pairs of identifiers.
  Converts identifiers using above function 'get_uniprot_id_mapping_pairs'
  but launches url requests in safe batches of 100 seqids at a time.

  from_type and to_type can be obtained from:
    http://www.uniprot.org/faq/28#mapping-faq-table
  """
  if cache_dir:
    seqids_txt = os.path.join(cache_dir, 'seqids.json')
    if os.path.isfile(seqids_txt):
        saved_seqids = json.load(open(seqids_txt))
        if saved_seqids != seqids:
            shutil.rmtree(cache_dir)

    if not os.path.isdir(cache_dir):
      os.makedirs(cache_dir)

    if not os.path.isfile(seqids_txt):
      json.dump(seqids, open(seqids_txt, 'w'))

  pairs = []
  i_seqid = 0
  if batch_size is None:
    batch_size = len(seqids)
  while i_seqid <= len(seqids):
    seqids_subset = seqids[i_seqid:i_seqid+batch_size]
    if cache_dir:
      subset_cache = os.path.join(cache_dir, 'mapping.%d.txt' % i_seqid)
    else:
      subset_cache = None
    subset_pairs = get_uniprot_id_mapping_pairs(
        from_type, to_type, seqids_subset, cache_fname=subset_cache)
    pairs.extend(subset_pairs)
    i_seqid += batch_size
  return pairs


def parse_isoforms(text):
  """
  Returns a dictionary of uniprot_acc and entries.

  Each entry contains:
    - var_seq: the sequence variations in the file
    - isoforms: the interpreted isoform sequences with
                the seqids
  """
  tag = None
  uniprot_data = {}
  var_seq = None
  in_isoform_section = False
  isoform_id = None
  for l in text.splitlines():
    test_tag = l[:5].strip()
    if test_tag and test_tag != tag:
      tag = test_tag
    line = l[5:].strip()
    words = line.split()
    if tag == "ID":
      uniprot_id = words[0]
      uniprot_data[uniprot_id] = {
        'var_seqs': [],
        'isoforms': {},
        'sequence': '',
      }
    if tag == "SQ":
      if words[0] != "SEQUENCE":
        uniprot_data[uniprot_id]['sequence'] += ''.join(words)
    if tag == 'FT':
      if var_seq is not None and l[5] != ' ':
        var_seq = None
      if line.startswith('VAR_SEQ'):
        var_seq = {
          'i': int(words[1]),
          'j': int(words[2]),
          'block': ''
        }
        uniprot_data[uniprot_id]['var_seqs'].append(var_seq)
      if var_seq is not None:
        var_seq['block'] += l[34:]
        if l.endswith('isoform'):
          # needed later to search isoform references
          var_seq['block'] += ' '
    if tag == 'CC':
      if words[0] == '-!-':
        if 'ALTERNATIVE' in words[1] and 'PRODUCTS' in words[2]:
          in_isoform_section = True
        else:
          in_isoform_section = False
      if in_isoform_section:
        if words[0].startswith('Name'):
          isoform_id = str(words[0][:-1].split('=')[1])
        for word in words:
          if word.startswith('IsoId='):
            seqid = word[:-1].split('=')[1]
            uniprot_data[uniprot_id]['isoforms'][isoform_id] = {
              'seqid': seqid
            }
            break
  for uniprot_id in uniprot_data:
    var_seqs = uniprot_data[uniprot_id]['var_seqs']
    isoforms = uniprot_data[uniprot_id]['isoforms']
    original_sequence = uniprot_data[uniprot_id]['sequence']
    for var_seq in var_seqs:
      block = var_seq['block']
      match = re.search(r'\(.*\)', block)
      isoform_tokens = match.group()[1:-1].split()
      isoform_ids = []
      for i in range(len(isoform_tokens)):
        if 'isoform' in isoform_tokens[i]:
          isoform_ids.append(str(isoform_tokens[i+1]))
      var_seq['isoform_ids'] = isoform_ids
      if block.startswith('Missing'):
        var_seq['deletion'] = True
      else:
        var_seq['deletion'] = False
        transition = block.split('(')[0]
        original, mutation = transition.split('->')
        var_seq['sequence'] = original.strip()
        assert len(var_seq['sequence']) == var_seq['j'] - var_seq['i'] + 1
        var_seq['mutation'] = mutation.strip()
    var_seqs.sort(key=lambda v:-v['i'])
    for isoform_id in isoforms:
      sequence = original_sequence
      for var_seq in uniprot_data[uniprot_id]['var_seqs']:
        if isoform_id in var_seq['isoform_ids']:
          i = var_seq['i']-1
          j = var_seq['j']
          if var_seq['deletion']:
            sequence = sequence[:i] + sequence[j:]
          else:
            sequence = sequence[:i] + var_seq['mutation'] + sequence[j:]
      isoforms[isoform_id]['sequence'] = sequence
  return uniprot_data


def parse_uniprot_txt_file(cache_txt):
  """
  Parses the text of metadata retrieved from uniprot.org.

  Only a few fields have been parsed, but this provides a
  template for the other fields.

  A single description is generated from joining alternative
  descriptions.

  Returns a dictionary with the main UNIPROT ACC as keys.
  """

  tag = None
  uniprot_id = None
  metadata_by_seqid = {}
  for l in cache_txt.splitlines():
    test_tag = l[:5].strip()
    if test_tag and test_tag != tag:
      tag = test_tag
    line = l[5:].strip()
    words = line.split()
    if tag == "ID":
      uniprot_id = words[0]
      is_reviewed = words[1].startswith('Reviewed')
      length = int(words[2])
      metadata_by_seqid[uniprot_id] = {
        'id': uniprot_id,
        'is_reviewed': is_reviewed,
        'length': length,
        'sequence': '',
        'accs': [],
      }
      entry = metadata_by_seqid[uniprot_id]
    if tag == "SQ":
      if words[0] != "SEQUENCE":
        entry['sequence'] += ''.join(words)
    if tag == "AC":
      accs = [w.replace(";", "") for w in words]
      entry['accs'].extend(accs)
    if tag == "DR":
      if 'PDB' in words[0]:
        if 'pdb' not in entry:
          entry['pdb'] = words[1][:-1]
        if 'pdbs' not in entry:
          entry['pdbs'] = []
        entry['pdbs'].append(words[1][:-1])
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
      if 'GO' in words[0]:
        if 'go' not in entry:
          entry['go'] = []
        entry['go'].append(' '.join(words[1:]))
      if 'Pfam' in words[0]:
        if 'pfam' not in entry:
          entry['pfam'] = []
        entry['pfam'].append(words[1][:-1])
    if tag == "GN":
      if 'gene' not in entry and len(words) > 0:
        pieces = words[0].split("=")
        if len(pieces) > 1 and 'name' in pieces[0].lower():
          entry['gene'] = pieces[1].replace(';', '').replace(',', '')
    if tag == "OS":
      if 'organism' not in entry:
        entry['organism'] = ""
      entry['organism'] += line
    if tag == "DE":
      if 'descriptions' not in entry:
        entry['descriptions'] = []
      entry['descriptions'].append(line)
    if tag == "CC":
      if 'comment' not in entry:
        entry['comment'] = ''
      entry['comment'] += line + '\n'

  for entry in metadata_by_seqid.values():
    descriptions = entry['descriptions']
    for i in reversed(range(len(descriptions))):
      description = descriptions[i]
      if 'Short' in description or 'Full' in description:
        j = description.find('=')
        descriptions[i] = description[j+1:].replace(';', '')
      else:
        del descriptions[i]
    entry['description'] = '; '.join(descriptions)

  return metadata_by_seqid


def parse_uniprot_metadata_with_seqids(seqids, cache_txt):
  """
  Returns a dictionary of metadata of given seqids, doing the requisite
  lookup of seqids in the ACCs and IDs, and more importantly handles
  isoform seqid's
  """
  metadata = parse_uniprot_txt_file(cache_txt)
  tmp = metadata.copy()
  for uniprot_id in metadata.keys():
    for seqid in metadata[uniprot_id]['accs'] + [metadata[uniprot_id]['id']]:
      tmp[seqid] = metadata[uniprot_id]
  metadata = tmp
  results = {}
  isoform_dict = parse_isoforms(cache_txt)
  for seqid in seqids:
    if seqid in metadata:
      results[seqid] = metadata[seqid]
    else:
      primary_seqid = seqid[:6]
      if primary_seqid in metadata:
        protein_metadata = metadata[primary_seqid]
        uniprot_id = protein_metadata['id']
        if uniprot_id in isoform_dict:
          isoforms = isoform_dict[uniprot_id]['isoforms'].values()
          for isoform in isoforms:
            if isoform['seqid'] == seqid:
              results[seqid] = deepcopy(protein_metadata)
              results[seqid]['accs'] = [seqid]
              results[seqid]['sequence'] = isoform['sequence']
  return results


def fetch_uniprot_metadata(seqids, cache_fname=None):
  """
  Returns a dictonary of the uniprot metadata (as parsed
  by parse_uniprot_txt_file) of the given seqids. The seqids
  must be valid uniprot identifiers.

  Now handles isoform versions of accession id's!
  """

  primary_seqids = [s[:6] for s in seqids]
  if cache_fname and os.path.isfile(cache_fname):
    logging("Loading cached metadata from " + cache_fname + "\n")
    cache_txt = open(cache_fname).read()
  else:
    logging("Fetching metadata for %d Uniprot IDs from http://uniprot.org ...\n" % len(primary_seqids))
    r = requests.post(
        'http://www.uniprot.org/batch/',
        files={'file':StringIO(' '.join(primary_seqids))},
        params={'format':'txt'})
    while 'Retry-After' in r.headers:
      t = int(r.headers['Retry-After'])
      logging('Waiting %d\n' % t)
      time.sleep(t)
      r = requests.get(r.url)
    cache_txt = r.text
    if cache_fname:
      open(cache_fname, 'w').write(r.text)
    if is_html(cache_txt):
      # Got HTML response -> error
      logging("Error in fetching metadata\n")
      return {}

  return parse_uniprot_metadata_with_seqids(seqids, cache_txt)


def batch_uniprot_metadata(seqids, cache_dir=None, batch_size=400):
  """
  Returns a dictonary of the uniprot metadata (as parsed
  by parse_uniprot_txt_file) of the given seqids. The seqids
  must be valid uniprot identifiers.

  Now handles isoform versions of accession id's!
  """

  unique_seqids = list(set(seqids))

  if cache_dir:
    seqids_txt = os.path.join(cache_dir, 'seqids.json')
    if os.path.isfile(seqids_txt):
        saved_seqids = json.load(open(seqids_txt))
        if saved_seqids != unique_seqids:
            shutil.rmtree(cache_dir)

    if not os.path.isdir(cache_dir):
      os.makedirs(cache_dir)

    if not os.path.isfile(seqids_txt):
      json.dump(unique_seqids, open(seqids_txt, 'w'))

  metadata = {}
  i_seqid = 0
  if batch_size is None:
    batch_size = len(unique_seqids)
  while i_seqid <= len(unique_seqids):
    seqids_subset = unique_seqids[i_seqid:i_seqid+batch_size]
    if cache_dir:
      subset_cache = os.path.join(cache_dir, 'metadata.%d.txt' % i_seqid)
    else:
      subset_cache = None
    metadata_subset = fetch_uniprot_metadata(
        seqids_subset, cache_fname=subset_cache)
    metadata.update(metadata_subset)
    i_seqid += batch_size
  return metadata


def is_text(seqid):
  if re.match('[A-Z,a-z,_]+$', seqid):
    return True
  return False


def is_uniprot(seqid):
  if re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]$', seqid):
    return True
  if re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]$', seqid):
    return True
  return False


def is_uniprot_variant(seqid):
  if is_uniprot(seqid[:6]):
    if len(seqid) == 6:
      return True
    variant = seqid[6:]
    if re.match('[-]\d+', variant):
      return True
  return False


def is_sgd(seqid):
  if re.match('Y[A-Z][L,R]\d\d\d[W|C]$', seqid):
    return True
  return False


def is_refseq(seqid):
  if re.match('[N,X,Y,Z][P,M]_\d+([.]\d+)?$', seqid):
    return True
  return False


def is_ensembl(seqid):
  if seqid.startswith("ENS"):
    return True
  return False


def is_maybe_uniprot_id(seqid):
  return '_' in seqid


def get_naked_seqid(seqid):
  if '|' not in seqid:
    return seqid
  pieces = seqid.split('|')
  if is_text(pieces[0]):
    return pieces[1]
  return seqid


assert is_refseq('NP_064308.1')
assert not is_refseq('NP_064308a1')
assert is_refseq('NP_064308')
assert is_sgd('YAL001C')
assert is_uniprot('A2AAA3')
assert not is_uniprot('A2AAA3-34')
assert is_uniprot_variant('A2AAA3-34')
assert is_uniprot('A2AAA3')
assert not is_uniprot_variant('A2AAA3-a')
assert not is_uniprot_variant('A2AAA3aaab')


def probe_id_type(entries, is_id_fn, name, uniprot_mapping_type, cache_fname):
  alternative_ids = []
  for entry in entries:
    if entry['id_type'] != '':
      continue
    if is_id_fn(entry['seqid']):
      entry['id_type'] = uniprot_mapping_type
      alternative_ids.append(entry['seqid'])
  if len(alternative_ids) == 0:
    return
  n_id = len(alternative_ids)
  pairs = batch_uniprot_id_mapping_pairs(
      uniprot_mapping_type, 'ACC', alternative_ids, cache_dir=cache_fname)
  alternative_to_uniprot = { p[0]:p[1] for p in pairs }
  for entry in entries:
    if entry['seqid'] in alternative_to_uniprot:
      entry['uniprot_acc'] = alternative_to_uniprot[entry['seqid']]


def get_metadata_with_some_seqid_conversions(seqids, cache_dir=None):
  logging("Looking up uniprot metadata for %d seqids\n" % len(seqids))
  entries = []
  for seqid in seqids:
    entries.append({
      'raw_seqid':seqid,
      'seqid':'',
      'id_type':'',
      'uniprot_acc':'',
      'metadata':''})

  # break up pieces when in form xx|xxxxx|xxx
  for entry in entries:
    entry['seqid'] = get_naked_seqid(entry['raw_seqid'])

  # convert a few types into uniprot_ids
  id_types = [
    (is_sgd, 'locustag', 'ENSEMBLGENOME_PRO_ID'),
    (is_refseq, 'refseqp', 'P_REFSEQ_AC'),
    (is_refseq, 'refseqnt', 'REFSEQ_NT_ID'),
    (is_ensembl, 'ensembl', 'ENSEMBL_ID'),
    (is_maybe_uniprot_id, 'uniprotid', 'ID')]
  for is_id_fn, name, uniprot_mapping_type in id_types:
    if cache_dir:
      seqid_cache_fname = os.path.join(cache_dir, name)
    else:
      seqid_cache_fname = None
    probe_id_type(entries, is_id_fn, name, uniprot_mapping_type, seqid_cache_fname)

  # delete the variant suffixes in some uniprot id's as probe_id_type
  # can't cope with uniprot variant for id mapping lookup (bad!)
  for entry in entries:
    if entry['id_type'] == '' and is_uniprot_variant(entry['seqid']):
      entry['seqid'] = entry['seqid'][:6]

  # map UNIPROT ID's to their current best entry
  if cache_dir:
    seqid_cache_fname = os.path.join(cache_dir, 'uniprotuniprot')
  else:
    seqid_cache_fname = None
  probe_id_type(entries, is_uniprot, 'UNIPROT-ACC', 'ACC+ID', seqid_cache_fname)

  uniprot_seqids = []
  for entry in entries:
    if 'uniprot_acc' in entry:
      uniprot_seqids.append(entry['uniprot_acc'])
    if is_uniprot_variant(entry['raw_seqid']):
      # put the isoform variants back up id name as
      # batch_uniprot_metadata can handle isoforms!
      uniprot_seqids.append(entry['raw_seqid'])
      entry['uniprot_acc'] = entry['raw_seqid']
  uniprot_dict = batch_uniprot_metadata(uniprot_seqids, cache_dir)

  result = {}
  for entry in entries:
    if entry['uniprot_acc'] == "":
      continue
    uniprot_acc = entry['uniprot_acc']
    if uniprot_acc not in uniprot_dict:
      continue
    result[entry['raw_seqid']] = uniprot_dict[uniprot_acc]

  return result


def get_filtered_uniprot_metadata(seqids, cache_txt):
  """
  Returns a dictionary of uniprot data, but filters
  seqids for uniprot identifiers by using a mapping call
  to uniprot first.
  """

  stripped_seqids = [s[:6] for s in seqids]
  pairs = batch_uniprot_id_mapping_pairs(
      'ACC+ID', 'ACC', stripped_seqids)
  uniprot_seqids = []
  for seqid1, seqid2 in pairs:
    if seqid1 in stripped_seqids and seqid1 not in uniprot_seqids:
      uniprot_seqids.append(seqid1)
  uniprot_dict = batch_uniprot_metadata(uniprot_seqids, cache_txt)
  for seqid in seqids:
    if seqid not in uniprot_seqids and seqid[:6] in uniprot_seqids:
      uniprot_dict[seqid] = uniprot_dict[seqid[:6]]
  return uniprot_dict


def sort_seqids_by_uniprot(seqids, uniprot_data):
  """
  Returns a sorted list of seqids such that longest proteins that are
  reviewed appear first in the list.
  """

  def cmp_longer_protein_is_first(seqid1, seqid2):
    return uniprot_data[seqid2]['length'] \
        - uniprot_data[seqid1]['length']

  def diff_list(orig_list, other_list):
    return [v for v in orig_list if v not in other_list]

  uniprot_seqids = filter(lambda s: s in uniprot_data, seqids)
  uniprot_seqids.sort(cmp=cmp_longer_protein_is_first)

  remainder_seqids = diff_list(seqids, uniprot_seqids)

  is_reviewed = lambda s: uniprot_data[s]['is_reviewed']
  reviewed_seqids = filter(is_reviewed, uniprot_seqids)
  reviewed_seqids.sort(cmp=cmp_longer_protein_is_first)

  unreviewed_seqids = diff_list(uniprot_seqids, reviewed_seqids)

  return reviewed_seqids + unreviewed_seqids + remainder_seqids


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
  if header.find("|") != -1 and '|' in header.split()[0]:
    tokens = header.split('|')
    # "gi|ginumber|gb|accession bla bla" becomes "gi|ginumber"
    seqid = "%s|%s" % (tokens[0], tokens[1].split()[0])
    name = seqid + ' ' + tokens[-1].strip()
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
  file to faciliate matching, but the keys in the returned
  structure corresponds to entries in 'seqids'.
  """
  live_seqid = None
  proteins = {}
  if seqid_fn is not None:
    original_seqid_map = { seqid_fn(s):s for s in seqids }
    seqids = original_seqid_map.keys()
  for i, line in enumerate(open(fasta_db)):
    if line.startswith(">"):
      fasta_seqid, description = \
          parse_fasta_header(line, seqid_fn)
      live_seqid = None
      for seqid in seqids:
        if fasta_seqid == seqid:
          live_seqid = fasta_seqid
          if seqid_fn:
            live_seqid = original_seqid_map[fasta_seqid]
          proteins[live_seqid] = {
            'sequence': "",
            'description': description,
          }
          break
    elif live_seqid:
      proteins[live_seqid]['sequence'] += line.strip()
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
    f.write(">" + seqid)
    if 'description' in proteins[seqid]:
      f.write(" " + proteins[seqid]['description'])
    f.write("\n")
    sequence = proteins[seqid]['sequence']
    for i in range(0, len(sequence), width):
      f.write(sequence[i:i+width] + "\n")
  f.close()
