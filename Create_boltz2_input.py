#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 19:45:07 2025

@author: joao
"""

import pandas as pd
import yaml
import os

df = pd.read_csv("Sequences_plus_smiles.csv")

def write_boltz_yaml(sequence, smiles, template_dir, output_name):
    """
    Escreve um arquivo YAML no formato esperado pelo Boltz, dado uma sequência de proteína e um SMILES.

    Args:
        sequence (str): sequência da proteína
        smiles (str): SMILES do ligante
        output_name (str): nome do arquivo de saída (padrão: 'input.yaml')

    Output:
        Cria um arquivo YAML em 'boltz_inputs/{output_name}'
    """
    # Criação do dicionário no formato correto
    data = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequence
                }
            },
            {
                "ligand": {
                    "id": "B",
                    "smiles": f"'{smiles}'"
                }
            }
        ],
        "templates" : [
            {
               "cif":  template_dir
            }
            ],
        "properties": [
            {
                "affinity": {
                    "binder": "B"
                }
            }
        ]
    }

    # Certifica-se que o diretório boltz_inputs existe
    os.makedirs("boltz_inputs", exist_ok=True)

    # Caminho do arquivo de saída
    output_path = os.path.join("boltz_inputs", output_name)

    # Escrita no arquivo YAML
    with open(output_path, "w") as f:
        yaml.dump(data, f, sort_keys=False)

    print(f"Arquivo YAML salvo em: {output_path}")
    
df.columns
template_path = "/home/joao/Documentos/Benchmark_boltz2/template/"
for _, row in df.iterrows():
    write_boltz_yaml(row["Sequência"], row["SMILES"],template_path, row["PDB ID"])
