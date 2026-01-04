#!/usr/bin/env python
"""Unit tests for the uniprot module."""

import unittest
import os
import tempfile
import textwrap
import uniprot


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TESTS_DIR, 'data')
EXAMPLES_DIR = os.path.join(os.path.dirname(TESTS_DIR), 'examples')


class TestSeqIdTypeDetection(unittest.TestCase):
    """Test sequence ID type detection functions."""

    def test_is_refseq_with_version(self):
        self.assertTrue(uniprot.is_refseq('NP_064308.1'))

    def test_is_refseq_without_version(self):
        self.assertTrue(uniprot.is_refseq('NP_064308'))

    def test_is_refseq_invalid(self):
        self.assertFalse(uniprot.is_refseq('NP_064308a1'))

    def test_is_sgd(self):
        self.assertTrue(uniprot.is_sgd('YAL001C'))

    def test_is_uniprot(self):
        self.assertTrue(uniprot.is_uniprot('A2AAA3'))

    def test_is_uniprot_invalid(self):
        self.assertFalse(uniprot.is_uniprot('A2AAA3-34'))

    def test_is_uniprot_variant(self):
        self.assertTrue(uniprot.is_uniprot_variant('A2AAA3-34'))

    def test_is_uniprot_variant_invalid_suffix_letters(self):
        self.assertFalse(uniprot.is_uniprot_variant('A2AAA3-a'))

    def test_is_uniprot_variant_invalid_suffix_mixed(self):
        self.assertFalse(uniprot.is_uniprot_variant('A2AAA3aaab'))


class TestFastaReading(unittest.TestCase):
    """Test FASTA file reading functionality."""

    def test_read_fasta_example(self):
        fasta_file = os.path.join(EXAMPLES_DIR, 'example.fasta')
        if os.path.exists(fasta_file):
            seqids, fastas = uniprot.read_fasta(fasta_file)
            self.assertIsInstance(seqids, list)
            self.assertIsInstance(fastas, dict)
            for seqid in seqids:
                self.assertIn(seqid, fastas)
                self.assertIn('sequence', fastas[seqid])
                self.assertIn('description', fastas[seqid])


class TestIsoformParsing(unittest.TestCase):
    """Test isoform parsing from UniProt metadata."""

    def setUp(self):
        self.isoform_dir = os.path.join(DATA_DIR, 'isoform')
        self.metadata_file = os.path.join(self.isoform_dir, 'Q91ZU6.txt')

    def test_parse_isoforms_structure(self):
        if os.path.exists(self.metadata_file):
            text = open(self.metadata_file).read()
            result = uniprot.parse_isoforms(text)
            self.assertIsInstance(result, dict)
            for uniprot_id, data in result.items():
                self.assertIn('var_seqs', data)
                self.assertIn('isoforms', data)
                self.assertIn('sequence', data)

    def test_isoform_sequences(self):
        if os.path.exists(self.metadata_file):
            seqids = ["Q91ZU6-{}".format(i) for i in [1, 2, 3, 4, 5, 6, 8]]
            text = open(self.metadata_file).read()
            results = uniprot.parse_uniprot_metadata_with_seqids(seqids, text)
            
            for seqid in seqids:
                fasta_file = os.path.join(self.isoform_dir, seqid + '.fasta')
                if os.path.exists(fasta_file):
                    read_seqids, fastas = uniprot.read_fasta(fasta_file)
                    test_sequence = list(fastas.values())[0]['sequence']
                    if seqid in results:
                        self.assertEqual(
                            test_sequence,
                            results[seqid]['sequence'],
                            f"Sequence mismatch for {seqid}"
                        )


class TestFastaWriting(unittest.TestCase):
    """Test FASTA file writing functionality."""

    def test_write_fasta_creates_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, 'test_output.fasta')
            proteins = {
                'seq1': {'sequence': 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV',
                         'description': 'Test protein 1'},
                'seq2': {'sequence': 'MTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV',
                         'description': 'Test protein 2'},
            }
            seqids = ['seq1', 'seq2']
            
            uniprot.write_fasta(output_file, proteins, seqids)
            
            self.assertTrue(os.path.exists(output_file))
            
            read_seqids, read_proteins = uniprot.read_fasta(output_file)
            self.assertEqual(read_seqids, seqids)
            for seqid in seqids:
                self.assertEqual(
                    read_proteins[seqid]['sequence'],
                    proteins[seqid]['sequence']
                )


class TestParseUniprotTxtFile(unittest.TestCase):
    """Test parsing UniProt metadata from text files."""

    def test_parse_isoform_metadata(self):
        metadata_file = os.path.join(DATA_DIR, 'isoform', 'Q91ZU6.txt')
        if os.path.exists(metadata_file):
            text = open(metadata_file).read()
            result = uniprot.parse_uniprot_txt_file(text)
            self.assertIsInstance(result, dict)
            
            for seqid, data in result.items():
                self.assertIn('id', data)
                self.assertIn('sequence', data)
                self.assertIn('is_reviewed', data)
                self.assertIn('length', data)


class TestParseFastaHeader(unittest.TestCase):
    """Test FASTA header parsing."""

    def test_parse_simple_header(self):
        header = "seq1 description here"
        seqid, name = uniprot.parse_fasta_header(header)
        self.assertEqual(seqid, 'seq1')
        self.assertIn('description', name)

    def test_parse_ncbi_header(self):
        header = ">gi|12345|gb|ABC123|protein description"
        seqid, name = uniprot.parse_fasta_header(header)
        self.assertEqual(seqid, 'gi|12345')

    def test_parse_header_with_angle_bracket(self):
        header = ">seq1 description"
        seqid, name = uniprot.parse_fasta_header(header)
        self.assertEqual(seqid, 'seq1')


class TestNakedSeqid(unittest.TestCase):
    """Test extracting naked seqid from formatted seqids."""

    def test_get_naked_seqid_with_pipes_text_first_gi(self):
        seqid = "gi|12345|gb|ABC123|description"
        naked = uniprot.get_naked_seqid(seqid)
        self.assertEqual(naked, '12345')

    def test_get_naked_seqid_with_pipes_text_first(self):
        seqid = "gb|ABC123|description"
        naked = uniprot.get_naked_seqid(seqid)
        self.assertEqual(naked, 'ABC123')

    def test_get_naked_seqid_without_pipes(self):
        seqid = "ABC123"
        naked = uniprot.get_naked_seqid(seqid)
        self.assertEqual(naked, seqid)


class TestEnsemblDetection(unittest.TestCase):
    """Test Ensembl ID detection."""

    def test_is_ensembl(self):
        self.assertTrue(uniprot.is_ensembl('ENSG00000196176'))
        self.assertFalse(uniprot.is_ensembl('Q91ZU6'))


class TestMaybeUniprotId(unittest.TestCase):
    """Test UniProt ID detection."""

    def test_is_maybe_uniprot_id(self):
        self.assertTrue(uniprot.is_maybe_uniprot_id('EFG_MYCA1'))
        self.assertFalse(uniprot.is_maybe_uniprot_id('A2AAA3'))


class TestUniProtMetadataFetching(unittest.TestCase):
    """Test UniProt metadata fetching."""

    def test_fetch_uniprot_metadata_with_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'metadata_cache.txt')
            
            metadata = textwrap.dedent("""\
                ID   TEST_PROT               Reviewed;         123 AA.
                AC   P12345;
                DE   RecName: Full=Test Protein;
                SQ   SEQUENCE                       10 AA;
                     MGTESTSEQ""")
            
            with open(cache_file, 'w') as f:
                f.write(metadata)
            
            result = uniprot.fetch_uniprot_metadata(['P12345'], cache_fname=cache_file)
            self.assertIsInstance(result, dict)


class TestSeqidtypeParsing(unittest.TestCase):
    """Test parsing of _SEQIDTYPE_SCRAPE."""

    def test_seqidtype_scrape_exists(self):
        self.assertTrue(hasattr(uniprot, '_SEQIDTYPE_SCRAPE'))
        self.assertGreater(len(uniprot._SEQIDTYPE_SCRAPE), 0)

    def test_get_seqidtype_id_types_returns_list(self):
        result = uniprot._get_seqidtype_id_types()
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)

    def test_seqidtype_id_types_are_strings(self):
        result = uniprot._get_seqidtype_id_types()
        for id_type in result:
            self.assertIsInstance(id_type, str)
            self.assertGreater(len(id_type), 0)


if __name__ == '__main__':
    unittest.main()

