#!/usr/bin/env python3
"""
Tests for the validation module.
"""

import sys
import os
import pytest

# Add parent directory to path so that we can import the scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.validation import validate_cat_order, parse_hmmscan_domtblout, extract_cat_module_sequences

def test_validate_cat_order_valid():
    """Test validate_cat_order with a valid C-A-T module."""
    domain_hits = [
        {'domain_type': 'C', 'start': 100, 'end': 400},
        {'domain_type': 'A', 'start': 420, 'end': 900},
        {'domain_type': 'T', 'start': 920, 'end': 1000}
    ]
    max_gap_allowed = 50
    
    is_valid, coords = validate_cat_order(domain_hits, max_gap_allowed)
    
    assert is_valid is True
    assert coords == (100, 1000)

def test_validate_cat_order_gap_too_large():
    """Test validate_cat_order with a gap that's too large."""
    domain_hits = [
        {'domain_type': 'C', 'start': 100, 'end': 400},
        {'domain_type': 'A', 'start': 500, 'end': 900},  # Gap of 100, which is > max_gap_allowed
        {'domain_type': 'T', 'start': 920, 'end': 1000}
    ]
    max_gap_allowed = 50
    
    is_valid, coords = validate_cat_order(domain_hits, max_gap_allowed)
    
    assert is_valid is False
    assert coords is None

def test_validate_cat_order_wrong_order():
    """Test validate_cat_order with domains in the wrong order."""
    domain_hits = [
        {'domain_type': 'A', 'start': 100, 'end': 400},
        {'domain_type': 'C', 'start': 420, 'end': 900},
        {'domain_type': 'T', 'start': 920, 'end': 1000}
    ]
    max_gap_allowed = 50
    
    is_valid, coords = validate_cat_order(domain_hits, max_gap_allowed)
    
    assert is_valid is False
    assert coords is None

def test_validate_cat_order_multiple_modules():
    """Test validate_cat_order with multiple C-A-T modules."""
    domain_hits = [
        {'domain_type': 'C', 'start': 100, 'end': 400},
        {'domain_type': 'A', 'start': 420, 'end': 900},
        {'domain_type': 'T', 'start': 920, 'end': 1000},
        {'domain_type': 'C', 'start': 1100, 'end': 1400},
        {'domain_type': 'A', 'start': 1420, 'end': 1900},
        {'domain_type': 'T', 'start': 1920, 'end': 2000}
    ]
    max_gap_allowed = 50
    
    is_valid, coords = validate_cat_order(domain_hits, max_gap_allowed)
    
    # Should return the first valid module
    assert is_valid is True
    assert coords == (100, 1000)

def test_validate_cat_order_overlapping_domains():
    """Test validate_cat_order with overlapping domains."""
    domain_hits = [
        {'domain_type': 'C', 'start': 100, 'end': 450},  # Overlaps with A
        {'domain_type': 'A', 'start': 420, 'end': 900},
        {'domain_type': 'T', 'start': 920, 'end': 1000}
    ]
    max_gap_allowed = 50
    
    is_valid, coords = validate_cat_order(domain_hits, max_gap_allowed)
    
    # This should still be valid since the domains are in the right order
    # and the gaps aren't too large (in fact, there's an overlap but not a gap)
    assert is_valid is True
    assert coords == (100, 1000)

def test_extract_cat_module_sequences():
    """Test extract_cat_module_sequences."""
    # Mock data
    hits_by_protein = {
        'protein1': [
            {'domain_type': 'C', 'start': 100, 'end': 400},
            {'domain_type': 'A', 'start': 420, 'end': 900},
            {'domain_type': 'T', 'start': 920, 'end': 1000}
        ]
    }
    sequences_by_id = {
        'protein1': 'A' * 100 + 'C' * 300 + 'G' * 20 + 'A' * 480 + 'T' * 20 + 'G' * 80  # 1000 aa
    }
    max_gap_allowed = 50
    
    module_sequences = extract_cat_module_sequences(hits_by_protein, sequences_by_id, max_gap_allowed)
    
    assert len(module_sequences) == 1
    assert 'protein1_CAT_Module_100_1000' in module_sequences
    assert len(module_sequences['protein1_CAT_Module_100_1000']) == 901  # 1000 - 100 + 1 