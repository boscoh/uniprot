from __future__ import print_function
import sys
import os
import uniprot
import json

scrape = """
UniProtKB AC/ID ACC+ID  from
UniProtKB AC  ACC to
UniProtKB ID  ID  to
UniParc UPARC both
UniRef50  NF50  both
UniRef90  NF90  both
UniRef100 NF100 both
Other sequence databases
EMBL/GenBank/DDBJ EMBL_ID both
EMBL/GenBank/DDBJ CDS EMBL  both
PIR PIR both
UniGene UNIGENE_ID  both
Entrez Gene (GeneID)  P_ENTREZGENEID  both
GI number*  P_GI  both
IPI P_IPI both
RefSeq Protein  P_REFSEQ_AC both
RefSeq Nucleotide REFSEQ_NT_ID  both
3D structure databases
PDB PDB_ID  both
DisProt DISPROT_ID  both
HSSP  HSSP_ID both
Protein-protein interaction databases
DIP DIP_ID  both
MINT  MINT_ID both
Protein family/group databases
Allergome ALLERGOME_ID  both
MEROPS  MEROPS_ID both
mycoCLAP  MYCOCLAP_ID both
PeroxiBase  PEROXIBASE_ID both
PptaseDB  PPTASEDB_ID both
REBASE  REBASE_ID both
TCDB  TCDB_ID both
PTM databases
PhosSite  PHOSSITE_ID both
Polymorphism databases
DMDM  DMDM_ID both
2D gel databases
Aarhus/Ghent-2DPAGE AARHUS_GHENT_2DPAGE_ID  both
World-2DPAGE  WORLD_2DPAGE_ID both
Protocols and materials databases
DNASU DNASU_ID  both
Genome annotation databases
Ensembl ENSEMBL_ID  both
Ensembl Protein ENSEMBL_PRO_ID  both
Ensembl Transcript  ENSEMBL_TRS_ID  both
Ensembl Genomes ENSEMBLGENOME_ID  both
Ensembl Genomes Protein ENSEMBLGENOME_PRO_ID  both
Ensembl Genomes Transcript  ENSEMBLGENOME_TRS_ID  both
GeneID  P_ENTREZGENEID  both
GenomeReviews GENOMEREVIEWS_ID  both
KEGG  KEGG_ID both
PATRIC  PATRIC_ID both
UCSC  UCSC_ID both
VectorBase  VECTORBASE_ID both
Organism-specific gene databases
AGD AGD_ID  both
ArachnoServer ARACHNOSERVER_ID  both
CGD CGD both
ConoServer  CONOSERVER_ID both
CYGD  CYGD_ID both
dictyBase DICTYBASE_ID  both
EchoBASE  ECHOBASE_ID both
EcoGene ECOGENE_ID  both
euHCVdb EUHCVDB_ID  both
EuPathDB  EUPATHDB_ID both
FlyBase FLYBASE_ID  both
GeneCards GENECARDS_ID  both
GeneFarm  GENEFARM_ID both
GenoList  GENOLIST_ID both
H-InvDB H_INVDB_ID  both
HGNC  HGNC_ID both
HPA HPA_ID  both
LegioList LEGIOLIST_ID  both
Leproma LEPROMA_ID  both
MaizeGDB  MAIZEGDB_ID both
MIM MIM_ID  both
MGI MGI_ID  both
neXtProt  NEXTPROT_ID both
Orphanet  ORPHANET_ID both
PharmGKB  PHARMGKB_ID both
PomBase POMBASE_ID  both
PseudoCAP PSEUDOCAP_ID  both
RGD RGD_ID  both
SGD SGD_ID  both
TAIR  TAIR_ID both
TubercuList TUBERCULIST_ID  both
WormBase  WORMBASE_ID both
WormBase Transcript WORMBASE_TRS_ID both
WormBase Protein  WORMBASE_PRO_ID both
Xenbase XENBASE_ID  both
ZFIN  ZFIN_ID both
Phylogenomic databases
eggNOG  EGGNOG_ID both
GeneTree  GENETREE_ID both
HOGENOM HOGENOM_ID  both
HOVERGEN  HOVERGEN_ID both
KO  KO_ID both
OMA OMA_ID  both
OrthoDB ORTHODB_ID  both
ProtClustDB PROTCLUSTDB_ID  both
Enzyme and pathway databases
BioCyc  BIOCYC_ID both
Reactome  REACTOME_ID both
UniPathWay  UNIPATHWAY_ID both
Gene expression databases
CleanEx CLEANEX_ID  both
GermOnline  GERMONLINE_ID both
Other
ChEMBL  CHEMBL_ID both
ChiTaRS CHITARS_ID  both
DrugBank  DRUGBANK_ID both
GenomeRNAi  GENOMERNAI_ID both
NextBio NEXTBIO_ID  both
"""
id_types = []
for line in scrape.splitlines():
  words = line.split()
  if words and words[-1] in ['both', 'to']:
    id_types.append(words[-2])


def analyze(seqid, cache_fname=None):
  if cache_fname is not None and os.path.isfile(cache_fname):
    cache = json.load(open(cache_fname))
  else:
    cache = { seqid: {} }
  print("===> Analyzing", seqid)
  good_types = []
  for from_type in id_types:
    if from_type not in cache[seqid]:
      pairs = uniprot.batch_uniprot_id_mapping_pairs(
        from_type, "ACC", [seqid])
      if pairs == []:
        cache[seqid][from_type] = None
      else:
        cache[seqid][from_type] = pairs[0][1]
      if cache_fname is not None:
        json.dump(cache, open(cache_fname, 'w'))
    if cache[seqid][from_type] is not None:
      good_types.append(from_type)
    print('{}:{} -> {}'.format(seqid, from_type, cache[seqid][from_type]))

  print(seqid, 'is compatible with:', ' '.join(good_types))


usage = """
seqidtype.py works out the type of seqid at uniprot.org by
brute-force matching a seqid against all seqid types

python seqidtype.py seqid1 seqid2 seqid3 ...

(Example: python seqidtype.py YOR261C)
"""


if __name__ == "__main__":
  if len(sys.argv) == 1:
    print(usage)
  else:
    for seqid in sys.argv[1:]:
      analyze(seqid, 'seqidtype.json')
