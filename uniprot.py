import os
import urllib
import urllib2
import pprint



# Doing internet lookups and converting is fraught with
# error, so a lot of caching to disk is recommended.
# Here's a couple of helper functions to cache Python
# dictionary data structures to disk, and back again


def read_dict(fname_dict):
  return eval(open(fname_dict).read())



def write_dict(a_dict, fname_dict, indent=2):
  pprint.pprint(
      a_dict, 
      indent=indent, 
      stream=open(fname_dict, 'w'))


def cache_dictionary_fn(fn, fname, is_overwrite=False):
  if not is_overwrite and os.path.isfile(fname):
    return read_dict(fname)
  else:
    results = fn()
    write_dict(results, fname)
    return results



def fetch_uniprot_id_mapping(fromtype, totype, seqids):
  """
  Converts a set of identifiers and returns a list of
  pairs of ids (original, new).

  fromtype and totype can be obtained from:
     http://www.uniprot.org/faq/28#id_mapping_examples
  """
  base = 'http://www.uniprot.org'
  tool = 'mapping'
  params = {
    'from': fromtype,
    'to': totype,
    'format': 'tab',
    'query': ' '.join(seqids),
  }
  data = urllib.urlencode(params)
  url = base+'/'+tool+'?'+data
  response = urllib2.urlopen(url)
  return response.read()



def batch_uniprot_id_mapping_pairs(
    fromtype, totype, seqids, n_batch=100):
  """
  Converts IDS using above function, but launches url
  requests in safe batches of 100 seqids at a time.
  Also filters results to produce a python list of
  matched pairs.
  """
  print "Fetching uniprot id mapping"
  txt = ""
  n_seqid = len(seqids)
  for i in range(0, n_seqid, n_batch):
    if n_seqid > n_batch:
      print "Looking up %d-%d seqid mappings" % (i, i+n_batch)
    txt += fetch_uniprot_id_mapping(fromtype, totype, seqids[i:i+n_batch])
  lines = filter(lambda l: 'from' not in l.lower(), txt.splitlines())
  return [l.split('\t')[:2] for l in lines]



def sequential_convert_ids_to_uniprot_id(seqids, fname_dict):
  """
  This is *the* function to convert arbitary seqid to a
  uniprot id. It only works for one seqid per url request, so
  is slow. But if you cannot transform your seqids so that
  batch mapping can work, it will be your only option. 

  Since this is slow, it caches the dictionary at every pass.
  """
  if not os.path.isfile(fname_dict):
    mapping = {}
  else:
    mapping = read_dict(fname_dict)
  for l in seqids:
    seqid = l.split()[0]
    if seqid not in mapping:
      result = batch_uniprot_id_mapping_pairs('ACC+ID', 'ACC', ['?' + seqid])
      if result:
        mapping[seqid] = result[0][0]
        write_dict(mapping, fname_dict)
        print seqid, "->", mapping[seqid]
      else:
        print seqid, '-> [null]'
  return mapping



def parse_uniprot_text(uniprot_fname):
  """
  Generates a dictionary of uniprot entries where the
  UNIPROT_ID is used as keys.

  Only a few fields have been parsed, but this provides a
  template for the other fields.

  Alternate descriptions are merged into one long string.
  """

  tag = None
  uniprot_id = None
  uniprot_data = {}
  for l in open(uniprot_fname):
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

  return uniprot_data



def fetch_uniprot_metadata_txt(seqids):
  uri = "http://uniprot.org/uniprot/"
  full_url = uri + "?"
  full_url += 'query='
  accs = [r'accession:' + s for s in seqids]
  full_url += r'+OR+'.join(accs)
  full_url += '&format=txt'
  return urllib2.urlopen(full_url).read()



def batch_uniprot_metadata(seqids, uniprot_fname, n_batch=100):
  """
  Returns a dictonary of the seqids that contains useful fields
  (as parsed by parse_uniprot_text), fetched from the
  www.uniprot.org. This uses batch processing in the url formation
  and seems only to work for less than 100 items at a time.
  """

  if os.path.isfile(uniprot_fname):
    print "Loading uniprot data from previous lookup", uniprot_fname
  else:
    print "Looking up %d seqids at http://www.uniprot.org" % len(seqids)
    try:
      f = open(uniprot_fname, 'w')
      for i in range(0, len(seqids), n_batch):
        print "Looking up %d-%d seqids" % (i, i+n_batch)
        txt = fetch_uniprot_metadata_txt(seqids[i:i+n_batch])
        f.write(txt)
      f.close()
    except:
      f.close()
      print "\nFailed connect to http://www.uniprot.org - limited metadata"
      return {}

  uniprot = parse_uniprot_text(uniprot_fname)

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



def parse_fasta_header(header, seqid_fn=None):
  """
  Parses a FASTA format header (with our without the initial '>') and returns a
  tuple of sequence id and sequence name/description.
  
  If NCBI SeqID format (gi|gi-number|gb|accession etc, is detected
  the first id in the list is used as the canonical id (see see
  http://www.ncbi.nlm.nih.gov/books/NBK21097/#A631 ).
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
  live_name = None
  proteins = {}
  if seqid_fn is not None:
    seqids = map(seqid_fn, seqids)
  for i, line in enumerate(open(fasta_db)):
    if line.startswith(">"):
      fasta_seqid, description = parse_fasta_header(line, seqid_fn)
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
  The main data-structure used in inmembrane is the proteins dictionary,
  which is initialized here from a source FASTA file. Also returns a list 
  of seqID's in the same order as that in the source FASTA file.
  """
  seqids = []
  seqid = None
  proteins = {}
  for line in open(fasta_db):
    if line.startswith(">"):
      seqid, name = parse_fasta_header(line, seqid_fn)
      seqids.append(seqid)
      proteins[seqid] = {
        'sequence':"",
        'name':name,
      }
      continue
    if seqid is not None:
      words = line.split()
      if words:
        proteins[seqid]['sequence'] += words[0]
  return seqids, proteins
  


