#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Created By  : Joao Paulo Ferreira
# Created Date: 04/10/2022
# version ='1.0'
# Module from Chloe: Blast methods
# ---------------------------------------------------------------------------
""" Este módulo será responsavel por fazer criação e buscas no banco de dados gerado blast."""
# ---------------------------------------------------------------------------
# Imports

import logging
import os

import docker
import pandas as pd


# ---------------------------------BLAST--------------------------------------

def create_blastdb(file):
    try:
        logging.getLogger("BLAST").info('Creating cluster database...')
        sintax = f"makeblastdb -in {file} -dbtype prot -input_type fasta"
        volume_file = f"{os.path.dirname(file)}:{os.path.dirname(file)}"
        logging.getLogger("BLAST").info(f"Blast db sintax: {sintax}")
        logging.getLogger("BLAST").info(f"Conteiner volume: {volume_file}")
        client = docker.from_env()
        client.containers.run("ncbi/blast:2.13.0", sintax, volumes=[volume_file])
        logging.getLogger("BLAST").info('Ending blast...')
    except Exception as erro:
        raise erro


def get_files_blast(path):
    try:
        return [file for file in os.listdir(path) if file.endswith(".blast.txt")]
    except Exception as erro:
        raise erro


class Blast:
    def __init__(self):
        self.table_columns = "\"6 qseqid qseq qstart qend sseqid sseq sstart send evalue bitscore score length pident "\
                             "nident mismatch gapopen gaps\" "

    def search_motif(self, path, query, db, name):
        try:
            """docker run -it -v $PWD:$PWD ncbi/blast blastp -query $PWD/teste.fa -db $PWD/dataset-prot.fa -out 
            $PWD/out.txt -sorthits 3 -num_threads 8 -outfmt "6 qseqid qseq qstart qend sseqid sseq sstart send evalue bitscore score length 
            pident nident mismatch gapopen gaps """
            logging.getLogger("BLAST_SEARCH").info("searching")
            sintax = f"blastp -sorthits 4 -num_threads 8 -gapopen 11 -qcov_hsp_perc 80 -query {query} -db {db} -out {path}/{name} -outfmt {self.table_columns}"
            volume = f"{path}:{path}"
            client = docker.from_env()
            logging.getLogger("BLASTP").info(sintax)
            client.containers.run("ncbi/blast:2.13.0", sintax, volumes=[volume])
        except Exception as erro:
            raise erro

    def parse_table_blast(self, path):
        try:
            columns = ['qseqid', 'qseq', 'qstart', 'qend', 'sseqid', 'sseq',
                       'sstart', 'send', 'evalue', 'bitscore', 'score', 'length',
                       'pident', 'nident', 'mismatch', 'gapopen', 'gaps', 'cluster']
            result = pd.DataFrame()
            for table in get_files_blast(path):
                df = pd.read_csv(os.path.join(path, table), sep='\t', header=None, index_col=None)
                df['cluster'] = table.split('.')[0]
                result = pd.concat([result, df], axis=0, ignore_index=True)
            result.columns = columns
            result.to_csv(f"{path}/RESULT.tsv", header=True, index=False, sep='\t')
        except Exception as erro:
            raise erro
