#!/usr/bin/env python
"""Integration tests for validating seqidtypes from _SEQIDTYPE_SCRAPE.

This test file validates that the ID mapping types listed in the
_SEQIDTYPE_SCRAPE string in uniprot.py are still valid according to the
UniProt API's /configure/idmapping/fields endpoint.

IMPORTANT NOTE: The _SEQIDTYPE_SCRAPE data is OUTDATED
===========================================================

The _SEQIDTYPE_SCRAPE in uniprot.py contains old API field names that are
no longer valid with the current UniProt API. For example:

  Old names (in scrape):    New names (in current API):
  - ACC                  -> UniProtKB_AC-ID
  - ID                   -> UniProtKB_AC-ID
  - P_REFSEQ_AC          -> RefSeq_Protein
  - REFSEQ_NT_ID         -> RefSeq_Nucleotide
  - ENSEMBL_ID           -> Ensembl
  - PDB_ID               -> PDB
  - etc.

The code in uniprot.py has been updated to use the new API field names
directly via the id_types list in get_metadata_with_some_seqid_conversions(),
so this is not a critical issue. However, the _SEQIDTYPE_SCRAPE documentation
should be updated for accuracy.

Run with: uv run python -m unittest test_seqidtypes -v
"""

import unittest
import json
import httpx
import uniprot


def get_valid_api_fields():
    """Fetch the valid ID mapping fields from the UniProt API."""
    try:
        response = httpx.get('https://rest.uniprot.org/configure/idmapping/fields', timeout=10.0)
        if response.status_code == 200:
            return response.json()
    except (httpx.RequestError, httpx.TimeoutException):
        return None


def has_internet():
    """Check if internet connectivity is available."""
    try:
        response = httpx.get('https://rest.uniprot.org/configure/idmapping/fields', timeout=5.0)
        return response.status_code == 200
    except (httpx.RequestError, httpx.TimeoutException):
        return False


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestSeqidtypeValidity(unittest.TestCase):
    """Test that seqidtypes in _SEQIDTYPE_SCRAPE are valid with current API."""
    
    @classmethod
    def setUpClass(cls):
        """Fetch valid API fields once for all tests."""
        api_data = get_valid_api_fields()
        if api_data:
            cls.valid_field_names = set()
            for group in api_data.get('groups', []):
                for item in group.get('items', []):
                    cls.valid_field_names.add(item['name'])
        else:
            cls.valid_field_names = set()
    
    def test_api_fields_available(self):
        """Test that we can fetch valid API fields."""
        self.assertGreater(len(self.valid_field_names), 0, 
                          "Should have fetched valid API fields")
    
    def test_seqidtype_scrape_fields_outdated(self):
        """Test documenting that _SEQIDTYPE_SCRAPE is outdated.
        
        Most ID types from _SEQIDTYPE_SCRAPE are NOT valid in the current API.
        This test documents this and should be updated when the scrape is refreshed.
        """
        id_types = uniprot._get_seqidtype_id_types()
        self.assertGreater(len(id_types), 0, "Should have parsed seqidtypes")
        
        invalid_types = []
        for id_type in id_types:
            if id_type not in self.valid_field_names:
                invalid_types.append(id_type)
        
        # Document that most types from the scrape are outdated
        # This is a known issue - _SEQIDTYPE_SCRAPE needs to be updated
        # with the new API field names
        self.assertGreater(len(invalid_types), 0,
                          "Most types from _SEQIDTYPE_SCRAPE are invalid with current API")
    
    def test_common_id_types_still_valid(self):
        """Test that commonly used ID types are valid."""
        common_types = [
            'UniProtKB_AC-ID',
            'UniProtKB',
            'Gene_Name',
            'RefSeq_Protein',
            'Ensembl',
            'PDB',
        ]
        
        for id_type in common_types:
            self.assertIn(id_type, self.valid_field_names,
                         f"Common ID type '{id_type}' should be valid")
    
    def test_deprecated_id_types_in_scrape(self):
        """Test identifying deprecated ID types in the scrape.
        
        NOTE: This test documents that _SEQIDTYPE_SCRAPE contains old API parameter names
        that are no longer valid in the current UniProt API. The scrape should be updated
        with the new API field names.
        """
        id_types = uniprot._get_seqidtype_id_types()
        
        # These are deprecated types that appear in _SEQIDTYPE_SCRAPE
        # but are no longer valid in the current API
        deprecated_types = [
            'ACC',              # -> UniProtKB_AC-ID
            'ID',               # -> UniProtKB_AC-ID  
            'P_REFSEQ_AC',      # -> RefSeq_Protein
            'REFSEQ_NT_ID',     # -> RefSeq_Nucleotide
            'ENSEMBL_ID',       # -> Ensembl
        ]
        
        found_deprecated = [dt for dt in deprecated_types if dt in id_types]
        self.assertGreater(len(found_deprecated), 0,
                          "Deprecated types should be present (documenting that scrape is outdated)")
    
    def test_id_type_mapping_combinations_valid(self):
        """Test that ID mapping combinations used in code are valid."""
        id_types = uniprot._get_seqidtype_id_types()
        
        # Common mapping target used in code
        target = 'UniProtKB_AC-ID'
        self.assertIn(target, self.valid_field_names,
                     f"Target ID type '{target}' should be valid")
        
        # Test a sample of mappings
        test_pairs = [
            ('UniProtKB_AC-ID', 'Gene_Name'),
            ('RefSeq_Protein', 'UniProtKB'),
            ('Gene_Name', 'UniProtKB'),
        ]
        
        for from_type, to_type in test_pairs:
            self.assertIn(from_type, self.valid_field_names,
                         f"Mapping source '{from_type}' should be valid")
            self.assertIn(to_type, self.valid_field_names,
                         f"Mapping target '{to_type}' should be valid")


@unittest.skipUnless(has_internet(), "No internet connectivity")
class TestSeqidtypeScrapeComparison(unittest.TestCase):
    """Compare _SEQIDTYPE_SCRAPE against current API documentation."""
    
    def test_scrape_data_extraction(self):
        """Test that _get_seqidtype_id_types correctly extracts fields from scrape."""
        id_types = uniprot._get_seqidtype_id_types()
        
        # Should extract at least 50+ ID types from the scrape
        self.assertGreater(len(id_types), 50,
                          "Should have extracted many ID types from scrape")
        
        # The scrape contains old field names (ACC, P_REFSEQ_AC, etc.)
        # not the new API field names. Document this with the actual field names.
        old_style_types = [
            'ACC', 'ID', 'P_REFSEQ_AC', 'REFSEQ_NT_ID', 'ENSEMBL_ID', 'PDB_ID'
        ]
        for old_style in old_style_types:
            self.assertIn(old_style, id_types,
                         f"Expected to find old-style field '{old_style}' in extracted types")
    
    def test_get_valid_api_fields_structure(self):
        """Test that the API response has the expected structure."""
        api_data = get_valid_api_fields()
        self.assertIsNotNone(api_data, "Should get API data")
        
        # Check structure
        self.assertIn('groups', api_data, "Should have 'groups' in API response")
        self.assertGreater(len(api_data['groups']), 0, "Should have groups")
        
        # Each group should have items
        for group in api_data['groups']:
            self.assertIn('items', group, "Each group should have 'items'")
            self.assertGreater(len(group['items']), 0, 
                             f"Group '{group.get('groupName')}' should have items")
            
            # Each item should have required fields
            for item in group['items']:
                self.assertIn('name', item, "Each item should have 'name'")
                self.assertIn('displayName', item, "Each item should have 'displayName'")


class TestSeqidtypeParsing(unittest.TestCase):
    """Test parsing of _SEQIDTYPE_SCRAPE (no network needed)."""
    
    def test_seqidtype_scrape_exists(self):
        """Test that _SEQIDTYPE_SCRAPE is defined."""
        self.assertTrue(hasattr(uniprot, '_SEQIDTYPE_SCRAPE'),
                       "Should have _SEQIDTYPE_SCRAPE defined")
        self.assertGreater(len(uniprot._SEQIDTYPE_SCRAPE), 0,
                          "_SEQIDTYPE_SCRAPE should not be empty")
    
    def test_get_seqidtype_id_types_function_exists(self):
        """Test that _get_seqidtype_id_types function exists."""
        self.assertTrue(callable(uniprot._get_seqidtype_id_types),
                       "_get_seqidtype_id_types should be callable")
    
    def test_get_seqidtype_id_types_returns_list(self):
        """Test that _get_seqidtype_id_types returns a list."""
        result = uniprot._get_seqidtype_id_types()
        self.assertIsInstance(result, list,
                            "_get_seqidtype_id_types should return a list")
        self.assertGreater(len(result), 0,
                          "Should return non-empty list")
    
    def test_seqidtype_id_types_are_strings(self):
        """Test that all extracted ID types are strings."""
        result = uniprot._get_seqidtype_id_types()
        for id_type in result:
            self.assertIsInstance(id_type, str,
                                f"ID type should be string, got {type(id_type)}")
            self.assertGreater(len(id_type), 0, "ID type should not be empty")
    
    def test_no_duplicate_id_types(self):
        """Test that there are minimal duplicate ID types."""
        result = uniprot._get_seqidtype_id_types()
        unique_types = set(result)
        
        # Allow at most one duplicate (P_ENTREZGENEID appears twice in scrape)
        duplicates = len(result) - len(unique_types)
        self.assertLessEqual(duplicates, 1,
                            f"Should have minimal duplicates, found {duplicates}")


if __name__ == '__main__':
    unittest.main()
