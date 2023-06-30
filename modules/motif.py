#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Created By  : Joao Paulo Ferreira
# Created Date: 15/05/2022
# version ='1.0'
# Module from Chloe: Motif Discovery
# ---------------------------------------------------------------------------
""" Este módulo será responsavel por fazer o descobrimento de motifs
em um dado conjunto de sequencias"""
# ---------------------------------------------------------------------------
# Imports

import logging
from Bio.motifs import meme
from Bio.SeqUtils import GC
import pandas as pd
import docker
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


# -----------------------------MOTIF DISCOVERY-------------------------------


def detect_ambiguous(seq):
    count = 0
    size = len(seq)
    for i in seq:
        if not i in ['A', 'C', 'T', 'G']:
            count = count + 1
    return f"{count}/{size}({round(((count * 100) / size), 2)}%)"


def write_fasta(output, seqs):
    try:
        handle = os.path.join(output, f"MOTIF_QUERY.fa")
        with open(handle, "a") as output_handle:
            SeqIO.write(seqs, output_handle, "fasta")

    except Exception as erro:
        raise erro


def df_to_fasta(df, output):
    try:
        records = list()
        for idx, row in df.iterrows():
            records.append(SeqRecord(Seq(row['MOTIF']),
                                     id=f"MOTIF-{idx + 1}-{row['CLUSTER']}",
                                     description=row['CLUSTER']))
        write_fasta(output, (record for record in records))
    except Exception as erro:
        raise erro


def parse_meme(file, output, len_dataset, base_path, id_cluster):
    try:
        logging.getLogger("PARSING_RESULTS").info("Parsing output of MEME")
        result = pd.DataFrame(columns=['MOTIF', 'STRAND', 'SITES',
                                       'PVALUE', 'GC', 'AMBIGUOUS_BASES', 'CLUSTER'])
        with open(file) as f:
            record = meme.read(f)

        for motif in record:
            for instance in motif.instances:
                if not result['MOTIF'].str.contains(instance.motif_name).any():
                    data = [result,
                            pd.DataFrame({
                                'MOTIF': [instance.motif_name],
                                'CLUSTER': id_cluster,
                                'STRAND': [instance.strand],
                                'SITES': [1],
                                'PVALUE': [instance.pvalue],
                                'GC': [f"{round((GC(instance.motif_name)), 2)}%"],
                                'AMBIGUOUS_BASES': [detect_ambiguous(instance.motif_name)]
                            })]
                    result = pd.concat(data, axis=0, ignore_index=True)
                else:
                    result.loc[result['MOTIF'] == instance.motif_name, 'SITES'] = \
                        result.loc[result['MOTIF'] == instance.motif_name, 'SITES'] + 1

        result['SITES'] = result['SITES'].map(
            lambda a: f"{a}/{len_dataset}({round(((a * 100) / len_dataset), 2)}%)")

        try:
            if not os.path.isfile(os.path.join(base_path, 'MOTIFS.tsv')):
                result.to_csv(os.path.join(base_path, 'MOTIFS.tsv'), sep='\t', index=False, header=True)
            else:
                result.to_csv(os.path.join(base_path, 'MOTIFS.tsv'), sep='\t', index=False, header=False, mode='a')

            logging.getLogger("PARSING_RESULTS").info(f"File saved: {os.path.join(base_path, 'MOTIFS.tsv')}")
        except Exception as erro:
            logging.getLogger("PARSING_RESULTS").error(f"Can't parse MEME output! {erro}")

        try:
            df_to_fasta(result, base_path)
        except Exception as erro:
            raise erro

    except Exception as erro:
        logging.getLogger("PARSE_MEME_RESULTS").error(erro)


class Motifs_discovery:
    def __init__(self, n_motifs, size_motif):
        self.n_motifs = n_motifs
        self.size_motif = size_motif

    def caller(self, file_input, output_path):
        try:
            logging.getLogger("MEME").info('Starting meme...')
            sintax = f"meme /home/meme/{file_input}  -p 8 -oc motif_disc -protein -nmotifs {self.n_motifs} -objfun " \
                     f"classic -w {self.size_motif} "
            volume = f"{str(output_path)}:/home/meme"
            logging.getLogger("MEME").info(f"Meme sintax: {sintax}")
            logging.getLogger('MEME').info(f"Conteiner volume: {volume}")
            client = docker.from_env()
            client.containers.run('memesuite/memesuite:5.5.0', sintax, volumes=[volume], cpu_count=8)
            logging.getLogger("MEME").info('Ending meme...')
        except Exception as erro:
            raise erro

    def main(self, file_input, output_path):
        try:
            self.caller(file_input, output_path)

        except Exception as erro:
            logging.getLogger('MEME').error(f"Error during execution:{erro}")
            raise erro
