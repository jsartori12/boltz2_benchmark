#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 17:08:59 2025

@author: joao
"""

import requests
from pathlib import Path
import pandas as pd
import os
from rdkit import Chem
from Bio import SeqIO
def download_pdb_and_fasta(pdb_code, base_output_dir="output_files"):
    """
    Baixa o arquivo PDB e o arquivo FASTA do RCSB e os salva em subpastas.

    Args:
        pdb_code (str): Código do PDB (ex: '1TUP').
        base_output_dir (str): Diretório base para as subpastas 'pdbs/' e 'fastas/'.
    """
    pdb_code = pdb_code.lower()
    base_path = Path(base_output_dir)
    pdb_dir = base_path / "pdbs"
    fasta_dir = base_path / "fastas"
    
    # Criar diretórios
    pdb_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    # URLs para download
    pdb_url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_code.upper()}/display"

    # Baixar PDB
    pdb_response = requests.get(pdb_url)
    if pdb_response.status_code == 200:
        pdb_path = pdb_dir / f"{pdb_code}.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_response.text)
        print(f"[✔] PDB salvo em: {pdb_path}")
    else:
        print(f"[✘] Erro ao baixar o arquivo PDB: {pdb_code.upper()}")

    # Baixar FASTA
    fasta_response = requests.get(fasta_url)
    if fasta_response.status_code == 200:
        fasta_path = fasta_dir / f"{pdb_code}.fasta"
        with open(fasta_path, "w") as f:
            f.write(fasta_response.text)
        print(f"[✔] FASTA salvo em: {fasta_path}")
    else:
        print(f"[✘] Erro ao baixar o arquivo FASTA: {pdb_code.upper()}")

download_pdb_and_fasta("3I7B")
df = pd.read_csv("TgCDPK1_PDBs.csv")
df.columns
for i in df.index:
    download_pdb_and_fasta(df.iloc[i,0])
    print(df.iloc[i,0])


def load_all_sdfs_from_directory(directory):
    all_mols = []
    ligans_name = []
    for filename in os.listdir(directory):
        if filename.endswith(".sdf"):
            filepath = os.path.join(directory, filename)
            suppl = Chem.SDMolSupplier(filepath)
            ligans_name.append(filename[:-4])
            for mol in suppl:
                if mol is not None:
                    all_mols.append(Chem.MolToSmiles(mol))
    data = {"Ligands": ligans_name,
            "SMILES": all_mols}
    data = pd.DataFrame(data)
    return data

directory = "caminho/para/seus/sdfs"
mol_list = load_all_sdfs_from_directory("/home/joao/Documentos/Benchmark_boltz2/output_files/ligands")

def read_fasta_sequence(fasta_path):
    """Lê a primeira sequência de um arquivo FASTA."""
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq)
    except Exception:
        return None
df_merged = df.merge(mol_list, how="left", left_on="Ligante principal", right_on="Ligands")
# Caminho dos FASTAs
fasta_dir = "output_files/fastas"
sequences = []

# Itera sobre cada linha do df_main
for _, row in df_merged.iterrows():
    pdb_id = row["PDB ID"].lower()
    
    # Lê a sequência FASTA
    fasta_path = os.path.join(fasta_dir, f"{pdb_id}.fasta")
    sequence = read_fasta_sequence(fasta_path)
    sequences.append(sequence)
    

# Adiciona ao novo DataFrame
df_combined = pd.DataFrame({
    "PDB ID": df["PDB ID"],
    "Sequência": sequences,
    "SMILES": df_merged["SMILES"]
})


df_combined = df_combined.merge(df[["PDB ID", "Cálcio"]], on="PDB ID", how="left")

df_combined.to_csv("Sequences_plus_smiles.csv", index = False)














