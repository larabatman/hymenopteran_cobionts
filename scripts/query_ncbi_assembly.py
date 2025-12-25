#!/usr/bin/env python3

"""
Query NCBI E-utilities to map an assembly accession (GCA_...) to basic metadata.
This script uses the core NCBI flow:
- esearch (assembly accession -> internal numeric UID)
- esummary (UID -> mtadata: organism, taxid, biosample, bioproject, ...)

Currently:
- Reads the first data row of Hymenopteran_genomes.csv
- Queries NCBI for that assembly
- Prints a single TSV line of results to stdout

This will later be extended to the full rows to write a full table
"""
from pathlib import Path
import csv
import json
import sys
import urllib.request
import urllib.parse

# Base URL for NCBI E-utilities endpoints
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def ncbi_esearch_assembly(gca):
    """
    This function converts the human-facing assembly accession (as GCA_964212075.1) into the NCBI internal numeric UID used by other endpoints.

    Parameters: 
    - gca: str, assembly accession 

    Returns: str | None, the first UID as a string if found, otherwise None. 
    """
    # Construct an NCBI search term targeting assembly accessions specifically
    term = f"{gca}[Assembly Accession]"
    # Query parameters for easeracg
    params = {
        "db": "assembly",  #which NCBI database to query
        "term": term,  #search expression
        "retmode": "json" #ask for json output which is easier to parse
    }
    # Build the URL: base + endpoint + URL-encoded parameters
    url = f"{NCBI_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params)
    # Perform HTTP GET, and parse json response
    with urllib.request.urlopen(url) as r:
        data = json.load(r)
    # esearch returns a list of matchin internal IDs under idlist
    ids = data["esearchresult"]["idlist"]
    # Return the first match if it exists and None otherwise
    return ids[0] if ids else None

def ncbi_esummary_assembly(uid):
    """
    This function retrieves assembly metadata from a numeric assembly UID. 
    Parameters:
    - uid: str, NCBI internal numeric UID returned by esearch

    Returns: dict, a dicitonary with key assembly metadata fields.
    """
    # Query parameters for esummary: uses the numeric UID!
    params = {
        "db": "assembly",
        "id": uid, 
        "retmode": "json"
    }
    # Build URL for esummary
    url = f"{NCBI_BASE}/esummary.fcgi?" + urllib.parse.urlencode(params)
    # Fetch and parse json
    with urllib.request.urlopen(url) as r: 
        data = json.load(r)
    # esummary result JSON: retrieve the list of returned UIDs as well as the actual record dictionary
    uid = data['result']['uids'][0] # This gets the UID stored as a string
    doc = data['result'][uid] # This gives all the metadata attached to that UID, as a dictionary from which we cna extract the fields from
    # Extract fields of interest
    return {
        "assembly_accession": doc.get("assemblyaccession"),
        "organism": doc.get("organism"),
        "taxid": doc.get("taxid"),
        "bioproject": doc.get("bioprojectid"),
        "biosample": doc.get("biosampleaccn"),
        "assembly_level": doc.get("assemblylevel")
    }

def main():
    """
    Main entry point to:
    - Read the first species row from Hympenopteran_genomes.csv
    - Query NCBI for assembly UID and metadata
    - Print results as TSV to stdout
    """
    # Path to script
    SCRIPT_DIR = Path(__file__).resolve().parent
    # Species table:
    infile = SCRIPT_DIR.parent / "species" / "Hymenopteran_genomes.csv"
    # Read CSV and grab the first data row
    with open(infile, newline = "") as f:
        reader = csv.DictReader(f)
        first_row = next(reader)
    # Extract fields from CSV header
    species = first_row['species']
    gca = first_row['accession']
    # Progress log to stderr so stdout is clean
    print(f"# Querying {species} ({gca})", file = sys.stderr)
    # Acession -> numeric UID
    uid = ncbi_esearch_assembly(gca)
    if uid is None:
        sys.exit(f"ERROR: no assembly UID for {gca}")

    # UID -> metadata
    meta = ncbi_esummary_assembly(uid)

    # Output header + one row as TSV
    print("species\tassembly_accession\ttaxid\tbioproject\tbiosample\tassembly_level")
    print(
        f"{species}\t"
        f"{meta['assembly_accession']}\t"
        f"{meta['taxid']}\t"
        f"{meta['bioproject']}\t"
        f"{meta['biosample']}\t"
        f"{meta['assembly_level']}"
    )

if __name__ == "__main__":
    main()