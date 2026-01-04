#!/usr/bin/env python
"""Integration tests for UniProt API.

These tests make actual HTTP calls to the UniProt API and should be run
separately from unit tests. They require network connectivity and may be slower.

Run with: uv run python -m unittest test_integration -v
"""

import unittest
import os
import tempfile
import time
import textwrap
import uniprot
import httpx


def has_internet():
    """Check if internet connectivity is available."""
    try:
        response = httpx.get('https://rest.uniprot.org/uniprotkb/search', 
                            params={'query': 'P12345', 'size': '1'},
                            timeout=5.0)
        return response.status_code in [200, 400]  # 400 means API is reachable
    except (httpx.RequestError, httpx.TimeoutException):
        return False


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestUniProtAPIIntegration(unittest.TestCase):
    """Integration tests with real UniProt API."""

    def test_id_mapping_acc_to_gene_name(self):
         """Test mapping UniProt accessions to gene names."""
         # Use a well-known protein - P69905 (Hemoglobin alpha)
         seqids = ['P69905']
         
         pairs = uniprot.get_uniprot_id_mapping_pairs(
             'UniProtKB_AC-ID', 'Gene_Name', seqids
         )
         
         self.assertGreater(len(pairs), 0, "Should find at least one mapping")
         self.assertEqual(pairs[0][0], 'P69905')
         # Gene name should be populated
         self.assertTrue(len(pairs[0][1]) > 0, "Gene name should not be empty")

    def test_id_mapping_refseq_to_uniprot(self):
         """Test mapping RefSeq protein to UniProt accession."""
         # Use a well-known RefSeq - NP_001005484 (human alpha-hemoglobin)
         seqids = ['NP_001005484']
         
         pairs = uniprot.get_uniprot_id_mapping_pairs(
             'RefSeq_Protein', 'UniProtKB', seqids
         )
         
         self.assertGreater(len(pairs), 0, "Should find RefSeq to UniProt mapping")
         self.assertEqual(pairs[0][0], 'NP_001005484')
         # Should map to a UniProt accession
         self.assertTrue(pairs[0][1].startswith('P') or pairs[0][1].startswith('Q'),
                        "Should map to UniProt accession")

    def test_batch_id_mapping(self):
        """Test batch ID mapping with multiple sequences."""
        # Use multiple well-known proteins
        seqids = ['P69905', 'P12345', 'Q9Y5K6']  # Hemoglobin, Calmodulin, etc.
        
        pairs = uniprot.batch_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', seqids, batch_size=10
        )
        
        self.assertGreater(len(pairs), 0, "Should find mappings")
        # Check that we got results for our input
        mapped_accs = [p[0] for p in pairs]
        self.assertIn('P69905', mapped_accs, "Should have mapping for P69905")

    def test_fetch_uniprot_metadata(self):
        """Test fetching metadata for a protein."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'metadata_cache.txt')
            
            # Fetch metadata for hemoglobin alpha
            seqids = ['P69905']
            metadata = uniprot.fetch_uniprot_metadata(seqids, cache_fname=cache_file)
            
            self.assertIn('P69905', metadata, "Should have metadata for P69905")
            entry = metadata['P69905']
            
            # Check required fields
            self.assertIn('sequence', entry)
            self.assertIn('length', entry)
            self.assertIn('is_reviewed', entry)
            
            # Hemoglobin should have sequence
            self.assertGreater(len(entry['sequence']), 0)
            
            # Verify cache was created
            self.assertTrue(os.path.exists(cache_file))

    def test_fetch_multiple_metadata(self):
        """Test fetching metadata for multiple proteins."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = os.path.join(tmpdir, 'metadata_cache')
            
            # Fetch metadata for multiple proteins
            seqids = ['P69905', 'P12345']  # Hemoglobin alpha, Calmodulin
            metadata = uniprot.batch_uniprot_metadata(seqids, cache_dir=cache_dir)
            
            # Should have entries for both
            self.assertGreater(len(metadata), 0)
            
            for seqid in seqids:
                if seqid in metadata:
                    entry = metadata[seqid]
                    self.assertIn('sequence', entry)
                    self.assertGreater(len(entry['sequence']), 0)

    def test_metadata_caching(self):
        """Test that metadata caching works correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'test_cache.txt')
            
            # First fetch - should create cache
            metadata1 = uniprot.fetch_uniprot_metadata(
                ['P69905'], cache_fname=cache_file
            )
            
            self.assertTrue(os.path.exists(cache_file), "Cache file should exist")
            
            # Second fetch - should read from cache
            metadata2 = uniprot.fetch_uniprot_metadata(
                ['P69905'], cache_fname=cache_file
            )
            
            # Results should be identical
            self.assertEqual(metadata1, metadata2)

    def test_parse_metadata_structure(self):
        """Test that parsed metadata has expected structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'metadata.txt')
            
            metadata = uniprot.fetch_uniprot_metadata(
                ['P69905'], cache_fname=cache_file
            )
            
            if 'P69905' in metadata:
                entry = metadata['P69905']
                
                # Check expected fields
                expected_fields = [
                    'id', 'is_reviewed', 'length', 'sequence', 'accs', 'description'
                ]
                for field in expected_fields:
                    self.assertIn(field, entry, 
                                 f"Metadata should contain '{field}' field")

    def test_isoform_metadata(self):
        """Test fetching protein with isoforms."""
        # P35557 has multiple isoforms
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'isoform_cache.txt')
            
            # Note: This test may take longer due to API processing
            seqids = ['P35557', 'P35557-1']  # Base + isoform
            
            metadata = uniprot.fetch_uniprot_metadata(
                seqids, cache_fname=cache_file
            )
            
            # Should find at least base protein
            self.assertGreater(len(metadata), 0)

    def test_seqidtype_analysis(self):
        """Test seqidtype analysis on real IDs."""
        # This is a slower test - use a single, simple ID
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = os.path.join(tmpdir, 'seqidtype_cache.json')
            
            # Test with a well-known SGD ID
            seqid = 'YOR261C'
            
            # This will make multiple API calls, so set a timeout
            # Note: This can take 30+ seconds as it tests all ID types
            import sys
            original_stdout = sys.stdout
            
            try:
                # Suppress logging output during test
                import io
                sys.stdout = io.StringIO()
                
                uniprot.seqidtype_analyze(seqid, cache_fname=cache_file)
                
                # Verify cache was created
                self.assertTrue(os.path.exists(cache_file),
                               "seqidtype cache should be created")
            finally:
                sys.stdout = original_stdout


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestAPIErrorHandling(unittest.TestCase):
    """Test error handling with real API."""

    def test_invalid_id_type_mapping(self):
        """Test mapping with invalid 'from' ID type."""
        # Using a fake ID type should fail gracefully
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'FAKE_ID_TYPE', 'Gene_Name', ['P69905']
        )
        
        # Should return empty list on error
        self.assertEqual(pairs, [])

    def test_empty_seqid_list(self):
        """Test mapping with empty sequence ID list."""
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', []
        )
        
        # Should handle gracefully
        self.assertIsInstance(pairs, list)

    def test_invalid_seqid_format(self):
        """Test mapping with invalid sequence ID."""
        # Using a sequence that doesn't exist should still return valid response
        pairs = uniprot.get_uniprot_id_mapping_pairs(
            'UniProtKB_AC-ID', 'Gene_Name', ['INVALID_ID_12345']
        )
        
        # Should return empty list (no matches) rather than error
        self.assertIsInstance(pairs, list)


if __name__ == '__main__':
    # Only run integration tests if explicitly requested
    import sys
    if 'test_integration' not in sys.argv[0]:
        print(textwrap.dedent("""\
            Integration tests require real API connectivity.
            Run with: uv run python -m unittest test_integration -v
            
            These tests will:
            - Make actual HTTP calls to UniProt API
            - Take several minutes to complete
            - Require stable internet connection
            """))
    
    unittest.main()
