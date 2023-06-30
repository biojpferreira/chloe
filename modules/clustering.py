#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Created By  : Joao Paulo Ferreira
# Created Date: 15/05/2022
# version ='1.0'
# Module from Chloe: Clustering
# ---------------------------------------------------------------------------
""" Este módulo será responsavel por realizar a clusterização do alinhamento."""
# ---------------------------------------------------------------------------
# Imports

import logging
from Bio import SeqIO
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Cluster import kmedoids
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


# ------------------------------CLUSTERING-----------------------------------

class Clustering:

    def __init__(self, alignment, output):
        self.output = output
        self.alignment = alignment  # arquivo de alinhamento
        self.samples = pd.DataFrame(columns=['ID', 'CLUSTER', 'SEQ'])  # dataframe que vai receber o ID do cluster
        self.sequences = list()  # Lista das sequencias a serem clusterizadas

    # Função responsável por receber uma lista das sequencias e realizar a clusterização
    def clustering(self):
        try:
            # TODO: Verificar a melhor forma de criar a matrix de distancia
            calculator = DistanceCalculator(model='identity')
            distance_matrix = calculator.get_distance(AlignIO.read(os.path.join(self.output, self.alignment), 'fasta'))
            # TODO: Adicionar parametro de numero de clusters
            clusterid, error, nfound = kmedoids(distance_matrix, nclusters=4, npass=1000)
            return clusterid
        except Exception as erro:
            raise erro

    # Função reponsável por receber a lista de clusteres e distrobuir nas respectivas amostras no DF
    def update_df(self, clusters):
        try:
            for i in range(len(clusters)):
                self.samples.iloc[i, self.samples.columns.get_loc('CLUSTER')] = clusters[i]
        except Exception as erro:
            raise erro

    # Função que cria diretório do cluster
    def create_dir(self, cluster_id):
        try:
            path_dir = os.path.join(self.output, f"cluster_{cluster_id}")
            if not os.path.exists(path_dir):
                os.makedirs(path_dir)
                return path_dir
        except Exception as erro:
            raise erro

    # Função responsável por escrever a sequencia no arquivo fasta
    def write_fasta(self, records):
        try:
            for cluster in records.keys():
                logging.getLogger('CLUSTERING').info(f"Size cluster {cluster}: {len(records[cluster])} sequences")
                path = self.create_dir(cluster)
                handle = os.path.join(path, f"cluster_{cluster}.fa")
                with open(handle, "w") as output_handle:
                    SeqIO.write(records[cluster], output_handle, "fasta")

        except Exception as erro:
            raise erro

    def create_fasta_cluster(self):
        try:
            records = {cluster: list() for cluster in self.samples["CLUSTER"].unique()}
            for idx, row in self.samples.iterrows():
                records[row['CLUSTER']].append(SeqRecord(Seq(row['SEQ'].replace('-', '')),
                                                         id=row['ID'],
                                                         description=f"cluster_{row['CLUSTER']}"))
            self.write_fasta(records)
        except Exception as erro:
            raise erro

    def check_size_cluster(self):
        try:
            for cluster in self.samples['CLUSTER'].unique():
                # TODO: atualizar valor minimo permitido por cluster
                if len(self.samples.loc[self.samples['CLUSTER'] == cluster]) <= 10:
                    logging.getLogger("CLUSTERING").error(f"Cluster {cluster} to small for analysis, please check "
                                                          f"parameters or sequences from input.")
                    exit()
        except Exception as erro:
            raise erro

    def main(self):
        try:
            for seq_record in SeqIO.parse(os.path.join(self.output, self.alignment), "fasta"):
                # Adiciona a sequencia na lista de sequencias
                self.sequences.append(str(seq_record.seq).upper())

                # Adiciona o ID da sequencia no DF de amostras
                self.samples = pd.concat([self.samples,
                                          pd.DataFrame(
                                              {"ID": [seq_record.id],
                                               "CLUSTER": [float('Nan')],
                                               "SEQ": [str(seq_record.seq).upper()]
                                               })],
                                         axis=0, ignore_index=True)

            self.update_df(self.clustering())

            # TODO: Verificar o tamanho do cluster e fechar por tamanho
            self.check_size_cluster()
            self.create_fasta_cluster()

        except Exception as erro:
            raise erro
