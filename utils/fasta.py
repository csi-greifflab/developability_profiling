"""
Work with fasta files
version 1.0

Author: Matteo Pariset

"""
import numpy as np

from Bio import SeqIO

import io
import requests

import sys

# WARNING: Using 'id' of SeqIO does NOT work for multiple chains (Biopython cuts it at the first whitespace)
def filter_chain(target_name):
    def _filter(chain):
        return np.any([x == target_name for x in chain.description.strip().split("|")[-1].split(" ")])
    return _filter

def get_sequence(seqs_db, pdb_id, chain_names):
    fasta_content = requests.get(f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id.lower()}/fasta").content.decode()
    if (pdb_id not in seqs_db) or not set(chain_names).issubset(seqs_db[pdb_id].keys()):
        for name in set(chain_names) - set(seqs_db.get(pdb_id, [])):
            candidates = list(filter(filter_chain(name), SeqIO.parse(io.StringIO(fasta_content), format="fasta")))
            if len(candidates) > 1:
                raise ValueError()
            elif len(candidates) == 0:
                print(f"Chain {name} not found", list(SeqIO.parse(io.StringIO(fasta_content), format="fasta")))
            
            seqs_db[pdb_id] = seqs_db.get(pdb_id, {})
            seqs_db[pdb_id].update({name: str(candidates[0].seq)})
        
        print(f"PDB: {pdb_id} \t {len(chain_names)} chains | inserted in seqs db")
    else:
        print(f"PDB: {pdb_id} \t {len(chain_names)} chains | skipped")

def dataframe_to_fasta(df):
    return "\n".join(df.apply(lambda x: f">{x[0]}\n{x[1]}", axis=1).values)

def save_dataframe_as_fasta(filepath, df):
    if not filepath.endswith(".fasta"):
        filepath += ".fasta"
    print(f"[INFO] Saving sequences in file {filepath}")
    with open(filepath, 'w') as f:
        f.write(dataframe_to_fasta(df))