#!/usr/bin/env python
"""Integration tests for UniProt API.

These tests make actual HTTP calls to the UniProt API and should be run
separately from unit tests. They require network connectivity and may be slower.

Run with: uv run python -m pytest tests/test_integration.py -v
"""

import unittest
import os
import tempfile
import textwrap
import uniprot
import httpx


def has_internet():
    """Check if internet connectivity is available."""
    try:
        response = httpx.get('https://rest.uniprot.org/uniprotkb/search', 
                            params={'query': 'P12345', 'size': '1'},
                            timeout=5.0)
        return response.status_code in [200, 400]
    except (httpx.RequestError, httpx.TimeoutException):
        return False


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestUniProtAPIIntegration(unittest.TestCase):
    """Integration tests with real UniProt API."""

    def test_id_mapping_acc_to_gene_name(self):
        seqids = ['P69905']
        
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', seqids
        )
        
        self.assertGreater(len(pairs), 0, "Should find at least one mapping")
        self.assertEqual(pairs[0][0], 'P69905')
        self.assertTrue(len(pairs[0][1]) > 0, "Gene name should not be empty")

    def test_id_mapping_refseq_to_uniprot(self):
        seqids = ['NP_001005484']
        
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'RefSeq_Protein', 'UniProtKB', seqids
        )
        
        self.assertGreater(len(pairs), 0, "Should find RefSeq to UniProt mapping")
        self.assertEqual(pairs[0][0], 'NP_001005484')
        self.assertTrue(pairs[0][1].startswith('P') or pairs[0][1].startswith('Q'),
                       "Should map to UniProt accession")

    def test_batch_id_mapping(self):
        seqids = ['P69905', 'P12345', 'Q9Y5K6']
        
        pairs = uniprot.batch_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', seqids, batch_size=10
        )
        
        self.assertGreater(len(pairs), 0, "Should find mappings")
        mapped_accs = [p[0] for p in pairs]
        self.assertIn('P69905', mapped_accs, "Should have mapping for P69905")

    def test_fetch_uniprot_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'metadata_cache.txt')
            
            seqids = ['P69905']
            metadata = uniprot.fetch_uniprot_metadata(seqids, cache_fname=cache_file)
            
            self.assertIn('P69905', metadata, "Should have metadata for P69905")
            entry = metadata['P69905']
            
            self.assertIn('sequence', entry)
            self.assertIn('length', entry)
            self.assertIn('is_reviewed', entry)
            
            self.assertGreater(len(entry['sequence']), 0)
            
            self.assertTrue(os.path.exists(cache_file))

    def test_fetch_multiple_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = os.path.join(tmpdir, 'metadata_cache')
            
            seqids = ['P69905', 'P12345']
            metadata = uniprot.batch_uniprot_metadata(seqids, cache_dir=cache_dir)
            
            self.assertGreater(len(metadata), 0)
            
            for seqid in seqids:
                if seqid in metadata:
                    entry = metadata[seqid]
                    self.assertIn('sequence', entry)
                    self.assertGreater(len(entry['sequence']), 0)

    def test_metadata_caching(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'test_cache.txt')
            
            metadata1 = uniprot.fetch_uniprot_metadata(
                ['P69905'], cache_fname=cache_file
            )
            
            self.assertTrue(os.path.exists(cache_file), "Cache file should exist")
            
            metadata2 = uniprot.fetch_uniprot_metadata(
                ['P69905'], cache_fname=cache_file
            )
            
            self.assertEqual(metadata1, metadata2)

    def test_parse_metadata_structure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'metadata.txt')
            
            metadata = uniprot.fetch_uniprot_metadata(
                ['P69905'], cache_fname=cache_file
            )
            
            if 'P69905' in metadata:
                entry = metadata['P69905']
                
                expected_fields = [
                    'id', 'is_reviewed', 'length', 'sequence', 'accs', 'description'
                ]
                for field in expected_fields:
                    self.assertIn(field, entry, 
                                 f"Metadata should contain '{field}' field")

    def test_isoform_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'isoform_cache.txt')
            
            seqids = ['P35557', 'P35557-1']
            
            metadata = uniprot.fetch_uniprot_metadata(
                seqids, cache_fname=cache_file
            )
            
            self.assertGreater(len(metadata), 0)

    def test_seqidtype_analysis(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'seqidtype_cache.json')
            
            seqid = 'YOR261C'
            
            import sys
            import io
            original_stdout = sys.stdout
            
            try:
                sys.stdout = io.StringIO()
                
                uniprot.seqidtype_analyze(seqid, cache_fname=cache_file)
                
                self.assertTrue(os.path.exists(cache_file),
                               "seqidtype cache should be created")
            finally:
                sys.stdout = original_stdout


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestAPIErrorHandling(unittest.TestCase):
    """Test error handling with real API."""

    def test_invalid_id_type_mapping(self):
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'FAKE_ID_TYPE', 'Gene_Name', ['P69905']
        )
        
        self.assertEqual(pairs, [])

    def test_empty_seqid_list(self):
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', []
        )
        
        self.assertIsInstance(pairs, list)

    def test_invalid_seqid_format(self):
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', ['INVALID_ID_12345']
        )
        
        self.assertIsInstance(pairs, list)


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestAPIFieldValidity(unittest.TestCase):
    """Validate that ID types used in the code are valid with current API."""

    @classmethod
    def setUpClass(cls):
        """Fetch valid API fields once for all tests."""
        try:
            response = httpx.get(
                'https://rest.uniprot.org/configure/idmapping/fields', 
                timeout=10.0
            )
            if response.status_code == 200:
                cls.valid_field_names = set()
                for group in response.json().get('groups', []):
                    for item in group.get('items', []):
                        cls.valid_field_names.add(item['name'])
            else:
                cls.valid_field_names = set()
        except (httpx.RequestError, httpx.TimeoutException):
            cls.valid_field_names = set()

    def test_common_id_types_valid(self):
        """Verify ID types used in the module are valid with current API."""
        common_types = [
            'UniProtKB_AC-ID',
            'UniProtKB',
            'Gene_Name',
            'RefSeq_Protein',
            'Ensembl',
        ]
        for id_type in common_types:
            self.assertIn(id_type, self.valid_field_names,
                         f"ID type '{id_type}' should be valid")


if __name__ == '__main__':
    import sys
    if 'test_integration' not in sys.argv[0]:
        print(textwrap.dedent("""\
            Integration tests require real API connectivity.
            Run with: uv run python -m pytest tests/test_integration.py -v
            
            These tests will:
            - Make actual HTTP calls to UniProt API
            - Take several minutes to complete
            - Require stable internet connection
            """))
    
    unittest.main()

