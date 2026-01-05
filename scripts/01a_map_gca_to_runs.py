#!/usr/bin/env python3

"""
Mapping from GCA to sequencing run with explicit evidence levels: taxid for species-level, bioproject for study-level and biosample for biological individual.

Input:
- A CSV file on STDIN with columns: 
    - species
    - accession (GCA_*)
Output:
- A TSV file on STDOUT
"""

# Imports and environment constraints
import csv #For CSV and TSV parsing
import json #For parsing JSON responses from NCBI and ENA APIs
import sys #For reading input from STDIN nad writing logs to STDERR
import time
import random #These are used for rate limiting and jitter, preventing APIs bans 
import socket #This was required for network-layer configuration
import urllib.parse
import urllib.request
import urllib.error #These libraries allow to build URLs, perform GET requests and catch any HTTP errors
# To work on the cluster, forcing IPv4 only was necessary as well as a 60s timeout
socket.setdefaulttimeout(60)
socket.has_ipv6 = False #Since many HPC nodes do not support IPv6, but urllib tries that first, it crashed with Errno 97

# Constants
# Base URLs for the two APIs to query
# NCBI containts the authoritative assembly metadata, and ENA the raw sequencing runs
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"
# Rate limit: wait 2s per request, which prevents HTTP 429: Too Many Requests errors
NCBI_SLEEP = 0.5

def ncbi_assembly_metadata(gca, retries = 5):
    """
    Resolve the GCA accession to assembly-level metadata: retreive what NCBI has explicitly declared.
    """
    # Retry loop since NCBI will occasionally fail, such that transient failurs do not affect the pipeline
    for attempt in range(1, retries + 1):
        try: 
            time.sleep(NCBI_SLEEP)
            # Step 1: esearch GCA to internal UID
            # Extract the NCBI internal numeric UIDs: GCAs are not used internally by NCBI, but everything works on numeric UIDs
            # Every database uses its own numeric primary keys, and esearch allows to retrieve it 
            # Build an HTTP GET request: create a params dict that urlencode serializes into a URL-encoded query string
            # The query is such as we search the assembly database looking for recordes where Assembly Accession == GCA_XXXX and return JSON
            params = {
                "db": "assembly",
                "term": f"{gca}[Assembly Accession]",
                "retmode": "json"
            }
            # Build the encoded GET URL
            url = f"{NCBI_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params)
            # Execute the request and parse the JSON in one step:
            with urllib.request.urlopen(url) as r:
                data = json.load(r) # data looks like {"esearchresult": {"count":1, "idlist": ["12345678"]}}
                ids = data.get("esearchresult", {}).get("idlist", []) # in "esearchresult", extract "idlist" and use an empty list if missing
                # If the GCA does not exist, stop
                if not ids:
                    return None
                # Use the first hit, as assemblies should be unique per GCA
                uid = ids[0] # This contains the internal asembly UID
                # Step 2: esummary from UID to the metadata
                time.sleep(NCBI_SLEEP)
                # Queyring the assembly summary record for assembly UID 
                # Not an search but a record lookup
                params = {
                    "db": "assembly",
                    "id": uid,
                    "retmode": "json"
                }
                url = f"{NCBI_BASE}/esummary.fcgi?" + urllib.parse.urlencode(params) 
                with urllib.request.urlopen(url) as r:
                    data = json.load(r) # data contains {"header": {...}, “result": {"uids": ["12345678"], "12345678": {"taxid": "36089", "bioproject": "PRJNA123456", "biosample": "SAMN98765432", ...}}} and so on
                    # The actual record is stored under a key named after the UID itself.
                doc = data["result"].get(uid, {}) # data["result"] is the dict holding all returned records. uid is the key specific to that assembly, and .get(uid, {}) retrieves that record or {} if missing
                # doc becomes doc = {"taxid": "362089", "bioproject": "PRHNA123456", "biosample": "SAMN98765432", ...} which is the assembly metadata object
                # Return the information of interest for the distance to raw sequencing data 
                return{
                    "taxid": doc.get("taxid"),
                    "bioproject": doc.get("bioproject"),
                    "biosample": doc.get("biosample")
                }
        # Error handling for HTTP failures explicitly
        except urllib.error.HTTPError as e:
            # Rate-limiting: back off
            if e.code == 429:
                wait = attempt * 5 + random.uniform(0, 3) # this is an exponential back off and jitter 
                print(f"# NCBI 429, sleeping {wait:.1f}s", file = sys.stderr)
                time.sleep(wait)
            else:
                 raise # any other HTTP error crashes
    return None

def ena_query(query):
    """
    Generic ENA query wrapper
    """
    # Query string passsed from main(): study_acession=PRJ... or sample_accession=SAM... or tax_id=XXXX
    params = {
        "result": "read_run", # Retrieveing runs
        "query": query,
        "fields": ",".join([ # Fields needed for classification: run UID, biosample, platform and strategy
            "run_accession",
            "sample_accession",
            "instrument_platform",
            "library_strategy"
        ]),
        "format": "json",
        "limit": 0 # this is ENA specific to return all results 
    }
    url = ENA_BASE + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout = 60) as r: 
        return json.load(r) # a list of dicts

def main():
    header = [
        "species",
        "assembly_accession",
        "taxid",
        "bioproject",
        "biosample",
        "run_mapping_level", # How we found the runs
        "has_pacbio",
        "has_hic",
        "has_rnaseq",
        "pacbio_runs",
        "hic_runs",
        "rnaseq_runs",
        "pacbio_biosamples",
        "hic_biosamples",
        "shared_pacbio_hic_biosample",
        "other_sequencing_methods" # What exists besides PacBio and HiC 
    ]
    # Read the raw CSV: each row becomes a dict keyed by the header names
    reader = csv.DictReader(sys.stdin)
    # Output TSV which is easier to parse downstream 
    writer = csv.writer(sys.stdout, delimiter = "\t", lineterminator = "\n")
    # Write the first row header 
    writer.writerow(header)

    # Per-row processing
    for row in reader:
        # One species/ assembly at a time, and taxid is always re-derived from NCBI
        species = row["species"]
        gca = row["accession"]

        print(f"# Processing {species} ({gca})", file = sys.stderr)
        # Anchor call 
        meta = ncbi_assembly_metadata(gca)
        # When the assembly is not found from the start, UNKNOWN row is emitted directly
        if not meta: 
            row = {col: "UNKNOWN" for col in header}
            row["species"] = species
            row["assembly_accession"] = gca
            writer.writerow(row[col] for col in header)
            continue
        taxid = meta.get("taxid")
        bioproject = meta.get("bioproject")
        biosample = meta.get("biosample")

        # Evidence hierarchy: keep track of how the runs are linked to the assembly 
        runs = []
        level = "NONE"

        try: 
            # The strongest link would be the same study as the assembly
            if bioproject:
                runs = ena_query(f"study_accession={bioproject}")
                level = "BIOPROJECT"
            # Fallback when study has no runs indexed in ENA
            if not runs and biosample:
                runs = ena_query(f"sample_accession={biosample}")
                level = "BIOSAMPLE"
            # Fallback when species-level is the only availability
            if not runs and taxid:
                runs = ena_query(f"tax_id={taxid}")
                level = "TAXID_FALLBACK"
        except Exception as e:
            print(f"# ENA error for {gca}: {e}", file = sys.stderr)

        # Run classification: keep a list of run accessions and sets to detect same-individual evidence
        pacbio_runs, hic_runs, rnaseq_runs = [], [], []
        pacbio_biosamples, hic_biosamples = set(), set()
        # Capture everything else 
        other_strategies = set()

        for r in runs:
            run = r.get("run_accession")
            bios = r.get("sample_accession")
            # Normalize metadata for consistency
            strat = (r.get("library_strategy") or "").upper()
            plat = (r.get("instrument_platform") or "").upper()

            if plat == "PACBIO_SMRT":
                pacbio_runs.append(run)
                if bios:
                    pacbio_biosamples.add(bios)
            elif strat in {"HI-C", "HIC"}:
                hic_runs.append(run)
                if bios:
                    hic_biosamples.add(bios)
            elif strat == "RNA-SEQ":
                rnaseq_runs.append(run)
            else:
                if strat:
                    other_strategies.add(strat)
    
        writer.writerow([
                species,
                gca,
                taxid or "UNKNOWN",
                bioproject or "UNKNOWN",
                biosample or "UNKNOWN",
                level,
                bool(pacbio_runs),
                bool(hic_runs),
                bool(rnaseq_runs),
                ",".join(sorted(set(pacbio_runs))),
                ",".join(sorted(set(hic_runs))),
                ",".join(sorted(set(rnaseq_runs))),
                ",".join(sorted(pacbio_biosamples)),
                ",".join(sorted(hic_biosamples)),
                bool(pacbio_biosamples & hic_biosamples),
                ",".join(sorted(other_strategies)) if other_strategies else ""
            ])
        time.sleep(2 + random.uniform(0, 1))

if __name__ == "__main__":
    main()