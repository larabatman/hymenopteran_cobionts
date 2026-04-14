 #!/usr/bin/env python3
"""
Annotate BioSample IDs with biological metadata (sex, tissue, developmental stage)

NCBI BioSample is used as the authoritive source. Since BioSamples (SAMMN, SAME) are not reliably efetchable, accession are converted to UID search, and the XML is parsed. 
"""
# Libraries
import csv # To preserve column order in TSVs
import json # To parse JSON responses from NCBI eserach
import sys # Access argv, stdin, stdout, stderr
import time 
import random # To throttle requests
import urllib.parse
import urllib.request # To build the HTTP GET requests
import xml.etree.ElementTree as ET # To parse the BioSample XML returned by efetch

# NCBI E-utilities base URL, where all NCBI API endpoints live
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI helper functions 
def esearch_biosample_uid(term):
    """
    Resolve a BioSample accession to an internal NCBI UID.

    Input: 
    - term (str): BioSample accession, as "SAMEA112975685"
    Ouput:
    - uid(str): internal numeric UID, as "34454967"
    or "" if nothing was found
    """
    # Parameters for esearch
    params = {
        "db": "biosample", # Target NCBI BioSample database
        "term": term, # Search term (accession)
        "retmode": "json" # JSON is easiest to parse
    }

    # Construct the full URL:
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)

    # Perform HTTP request
    with urllib.request.urlopen(url, timeout = 60) as r:
        data = json.load(r)
        # data structure: data["esearchresult"]["idlist"] = list of matching UIDs in string
    ids = data.get("esearchresult", {}).get("idlist", [])

    # By convention, take the first UID if present 
    return ids[0] if ids else ""

def efetch_biosample_xml(uid):
    """
    Fetch BioSample XML given a UID. 
    efetch(db=biosample, id=SAMEA...) often returns empry <BioSampleSet/>, but efetch(db=biosample, id=UID) always returns the record

    Input: 
    - uid (str): internal NCBI BioSample UID
    - root (Element): root of parsed XML tree
    """
    params = {
        "db": "biosample", # BioSample database
        "id": uid, # Internal UID
        "retmode": "xml" # BioSample metadata is XML-only
    }
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params)

    with urllib.request.urlopen(url, timeout = 60) as r:
        tree = ET.parse(r)
        return tree.getroot() # Typically <BioSampleSet> or <BioSample>
    
def ncbi_biosample_metadata(sample_acc, visited = None, retries = 3):
    """
    Resolve biological metadata for a BioSample accession. 
    This function allws to resolve accession to UID, fetches XML and extracts sex/ tissue/ dev_stage. It follows parent samples while avoiding infinite recursions.

    Input:
    - sample_acc (str): SAMEA... or SAMN...
    - visited (set): accessions already visited for loop protection
    - retries (int): how many times to retry on transient failure
    Output:
    - dict with keys sex, tissue, dev_stage, description 
    or {} if nothing could be resolved
    """
    # Initialize visited set on first call
    if visited is None:
        visited = set()
    # Loop protection: if already visited accession, stop
    if sample_acc in visited:
        return {}
    visited.add(sample_acc)

    # Step 1: accession to UID
    uid = esearch_biosample_uid(sample_acc)
    if not uid:
        # If NCBI cannot resolve it, done
        return {}
    
    # Retry loop for robustness
    for _ in range(retries):
        try: 
            # Step 2: fetch XML by UID
            root = efetch_biosample_xml(uid)

            meta = {} # collected biological metadata
            parent = None # parent BioSample, if exist

            # First pass: detect parent BioSample which is often encoded as an Attribute called "sample same as"
            for attr in root.findall(".//{*}Attribute"):
                if attr.attrib.get("attribute_name", "").lower() == "sample same as":
                    parent = (attr.text or "").strip()
            # Second pass: extract biological attributes with {*} that majes it namespace-safe
            for attr in root.findall(".//{*}Attribute"):
                key = attr.attrib.get("attribute_name", "").lower()
                key = key.replace("_", " ").strip()
                value = (attr.text or "").strip()

                if not value:
                    continue

                if key == "sex":
                    meta["sex"] = value 
                
                elif key in {
                    "tissue",
                    "tissue type",
                    "organism part", 
                    "organism_part"
                }:
                    meta["tissue"] = value
                
                elif key in {
                    "life stage",
                    "dev stage",
                    "developmental stage",
                    "lifestage"
                }:
                    meta["dev_stage"] = value
                
            # Description/ title 
            title = root.find(".//{*}Title")
            if title is not None and title.text:
                meta["description"] = title.text.strip()
            
            # If anything biological is found, return it
            if meta:
                return meta
            # Otherwise follow the parent BioSample, when present
            if parent:
                return ncbi_biosample_metadata(parent, visited = visited)
            
            # Nothing found at all:
            return {}
        except Exception:
            # Transient failure: sleep and retry
            time.sleep(2 + random.uniform(0, 1))
    return {}

# TSV helper
def split_samples(value):
    """
    Split comma-separated BioSample accessions from TSV fields
    
    Input:
    - value (str): as "SAMEA1,SAMEA2"
    Output:
    - list[str]: ["SAMEA1", "SAMEA2"]
    """
    return [v.strip() for v in value.split(",") if v.strip()]

def pick(*values):
    """
    Return the first non-empty value. 

    Used to safely handle missing metadata fields.
    """
    for v in values:
        if v:
            return v
    return ""

# Main logic
def main():
    # Input TSV: either file passed as argument or STDIN
    infile = sys.argv[1]
    fin = open(infile, newline = "")

    # Read TSV rows as dicts keyed by column name
    reader = csv.DictReader(fin, delimiter = "\t")
    # Writer preserves tab separation and column order
    writer = csv.writer(sys.stdout, delimiter = "\t", lineterminator = "\n")
    # Output header: original columns + annotation columns
    writer.writerow(reader.fieldnames + [
        "sample_accession",
        "sex",
        "tissue",
        "dev_stage",
        "description"
    ])

    # Cache to avoide querying the same BioSample multiple times:
    cache = {}

    for row in reader:
        species = row.get("species", "UNKNOWN")
        print(f"# Processing species: {species}", file = sys.stderr)

        # Collect BioSample accessions from PacBio and Hi-C columns
        samples = set(
            split_samples(row.get("pacbio_biosamples", ""))+
            split_samples(row.get("hic_biosamples", ""))
        )

        for s in samples:
            print(f"#   Resolving BioSample: {s}", file = sys.stderr)

            # Query NCBI only once per BioSample
            if s not in cache:
                cache[s] = ncbi_biosample_metadata(s)
                time.sleep(1 + random.uniform(0, 1))
            
            meta = cache[s]

            # Emit enriched row
            writer.writerow([
                *[row[f] for f in reader.fieldnames], # Keep all the rows from input
                s,
                pick(meta.get("sex")) or "UNKNOWN",
                pick(meta.get("tissue")) or "UNKNOWN",
                meta.get("dev_stage", ""),
                meta.get("description", "")
            ])
        # If no BioSamples at all, emit the row with empty annotation
        if not samples:
            writer.writerow([
                *[row[f] for f in reader.fieldnames],
                "",          # sample_accession
                "UNKNOWN",   # sex
                "UNKNOWN",   # tissue
                "",          # dev_stage
                ""           # description
            ])
    if infile != "-":
        fin.close()
if __name__ == "__main__":
    main()