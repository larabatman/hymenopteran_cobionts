#!/usr/bin/env python3

"""
This script aims at determining what raw sequencing data exists for a species associated with a GCA. 

- GCA is used to resolve taxid of the assembly
- Raw data is discovered independently, via taxid (ENA)
- The coherence is then assessed via BioSample overlap.

This way, all unreliable assembly to SRA linkages are avoided. 
"""

from pathlib import Path # to build filesystem paths
import csv # to read the species genomes file
import json # to parse tthe NCBI/ ENA responses
import sys 
import urllib.request # to perform HTTP GET requests
import urllib.parse # to encode the URL query parameters
import time # to add some waiting time in between requests
import random
import urllib.error
# Force IPv4 networking because HPC environment does not support IPv6, and sometimes urllib fails because of that.
import socket
socket.setdefaulttimeout(60)
socket.has_ipv6 = False

# Constant used for NCBI and ENA querying:
# These are the base URL roots for the APIs, onto which to append endpoint and encoded parameters
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"

# Helper functions for NCBI-related requests
def ncbi_assembly_metadata(gca):
    """
    This function aims at resolving a GCA accession to:
    - taxid for the organism identity
    - BioProject for the conceptual study linkage
    The assembly itself is not used to enumerate reads: it is used as a reference. 
    
    Input: gca, str
    Output: assembly bioproject, str 
    """
    time.sleep(0.34) # 3 requests/s, to allow NCBI to breath
    # GCA towards internal assembly UID: 
    # easearch
    params = {
        "db": "assembly", #search in NCBI Assembly database
        "term": f"{gca}[Assembly Accession]", #query must match accesstion
        "retmode": "json"
    }
    # Use urllib.parse.urlencode(params) to turn the dict into db=assembly&term=...=&retmode=json, and url becomes full GET URL
    url = f"{NCBI_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params) 
   
    with urllib.request.urlopen(url) as r: #r: HTTP response file-like object
        data = json.load(r) # parse entire response into a dict: data
        # data structure: data["esearchresult"]["idlist"] is a list of internal NCBI numeric UIDs, as strings

    idlist = data["esearchresult"]["idlist"]
    if not idlist:
        return None

    uid = idlist[0] # this is the string numeric ID, not the GCA but NCBI's internal assembly record ID
    
    # UID towards metadata
    # esummary
    params = {
        "db": "assembly",
        "id": uid,
        "retmode": "json"
    }
    # query the metadata about that assembly UID
    url = f"{NCBI_BASE}/esummary.fcgi?" + urllib.parse.urlencode(params) 
    
    with urllib.request.urlopen(url) as r:
        data = json.load(r)
        # data["result"] is a dict mapping uid to record dict, and contains many fields: doc["taxid"], doc["biosampleaccn"] and so on
    doc = data["result"][uid]

    return{
        "taxid": doc.get("taxid")
    }

# Helper functions for ENA-related requests:
def ena_runs_from_taxid(taxid, max_retries = 3, base_sleep = 5):
    """
    This function aims to retrieve all sequencing runs for a species via ENA. 
    This helps define what data exists at all.

    Input: taxid as str or int
    Output: list[dict], where each dict is a run record
    """
    params = {
    "result": "read_run",
    "query": f"tax_id={taxid}",
    "fields": ",".join([
        "run_accession",
        "sample_accession",
        "instrument_platform",
        "library_strategy"
    ]),
    "format": "json",
    "limit": 0
    }

    # Build the GET URL for ENA
    url = ENA_BASE + "?" + urllib.parse.urlencode(params)
    
    for attempt in range(1, max_retries + 1):
        try: 
            with urllib.request.urlopen(url) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            # ENA wants to slow down, or request too big
            print(
                f"# ENA HTTP {e.code} for taxid {taxid}"
                f"(attempt {attempt}/{max_retries})",
                file = sys.stderr
            )
        # Add exponential backoff and jitter
        sleep_time = base_sleep * attempt + random.uniform(0, 2)
        time.sleep(sleep_time)
    # Give up otherwise
    return []

# Helper
def uniq(values):
    """
    This function returns sorted, unique, non-empty values.
    """
    return sorted({v for v in values if v})
            
# Main 
def main():
    # Read one species row for now
    script_dir = Path(__file__).resolve().parent # this makes an absolute path to the file containing this script
    infile = script_dir.parent / "species" / "Hymenopteran_genomes.csv" #this moves up one directory, and appends the path to the csv
    
    # Use the writer functions from csv to avoid discrepancies in column counts
    writer = csv.writer(
        sys.stdout, 
        delimiter = "\t",
        lineterminator = "\n"
    )
    # Print header once 
    writer.writerow([
            "species",
            "assembly_accession",
            "taxid",
            "has_pacbio",
            "has_hic",
            "has_rnaseq",
            "pacbio_runs",
            "hic_runs",
            "rnaseq_runs",
            "pacbio_biosamples",
            "hic_biosamples",
            "shared_pacbio_hic_biosample",
    ])

    with open(infile, newline="") as f: #open the csv
        reader = csv.DictReader(f) #yield each row as dict using header keys
        # row is dict like {"species": "Abia_cadens", "accession": "GCA_...", "taxid":"362089", ...}
    
        for row in reader: 
            species = row["species"] # extract species
            gca = row["accession"] # extract GCA

            print(f"# Querying {species} ({gca})", file = sys.stderr) #print to stderr to keep stdout TSV clean

            # 1) GCA to taxid and Bioproject
            taxid = None
            pacbio, hic, rnaseq = [], [], []
            try: 
                meta = ncbi_assembly_metadata(gca) # Call NCBI, and retrieves taxid and assembly project
                if meta is None:
                    raise ValueError("assembly not found")
                taxid = meta.get("taxid")
                if not taxid:
                        print(f"# No taxid for {gca}", file = sys.stderr)
                        taxid = None
            except Exception as e:
                print(f"# NCBI failed for {gca}: {e}", file = sys.stderr)
            
            # 2) taxid to all ENA runs
            runs = []
            if taxid:
                try:
                    runs = ena_runs_from_taxid(taxid) # Calls ENA, and builds the list of dicts for each run
                except Exception as e:
                    print(f"# ENA failed for {gca}: {e}", file = sys.stderr)
            
            # 3) classify and track linkage
            for r in runs: # r is a dict from ENA
                entry = {
                    "run": r["run_accession"], # ERR...
                    "biosample": r.get("sample_accession") #SAMEA...
                }
                platform = r.get("instrument_platform", "")
                strategy = r.get("library_strategy", "")
                if platform == "PACBIO_SMRT":
                    pacbio.append(entry)
                if strategy.upper() in {"HI-C", "HIC"}: # Note: case-sensitivity, if HI-C will be missed 
                    hic.append(entry)
                if strategy.upper() == "RNA-SEQ":
                    rnaseq.append(entry)
            pacbio_runs = uniq(e["run"] for e in pacbio)    
            hic_runs = uniq(e["run"] for e in hic)  
            rnaseq_runs = uniq(e["run"] for e in rnaseq)  

            pacbio_samples = uniq(e["biosample"] for e in pacbio)
            hic_samples = uniq(e["biosample"] for e in hic)

            shared_sample = bool(set(pacbio_samples) & set(hic_samples))

            writer.writerow([
                    species,
                    gca,
                    taxid or "",
                    bool(pacbio_runs),
                    bool(hic_runs),
                    bool(rnaseq_runs),
                    ",".join(pacbio_runs),
                    ",".join(hic_runs),
                    ",".join(rnaseq_runs),
                    ",".join(pacbio_samples),
                    ",".join(hic_samples),
                    shared_sample
            ])
            # Being polite and waiting before next species:
            time.sleep(3 + random.uniform(0, 2))

if __name__ == "__main__":
    main()