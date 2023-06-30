#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Created By  : Joao Paulo Ferreira
# Created Date: 15/05/2022
# version ='1.0'
# ---------------------------------------------------------------------------
""" Clean, align and clustering sequences for search for conserved regions inside."""
# ---------------------------------------------------------------------------
# Imports

import argparse
import os
import sys
import datetime
import string
import secrets
import logging
import concurrent.futures
import stat
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Imports internal modules
from modules import sequence_cleaner
from modules import alignment
from modules import clustering
from modules import motif
from modules import blast

# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', required=True, help="Path for fasta file.")
parser.add_argument('--output', '-o', required=True, help="Path job output.")
parser.add_argument('--min_lenght', type=float, default=0, required=False,
                    help="""the user defines the minimum length. DEFAULT = 0, it means you don’t have 
                    to care about the minimum length""")

parser.add_argument('--max_lenght', type=float, default=0, required=False,
                    help="""the user defines the maximum length. DEFAULT = 0, it means you don’t have 
                    to care about the maximum length""")

parser.add_argument("--percent_n", required=False, type=float, default=100, help=""" The user defines the percent of N is allowed. 
                    DEFAULT = 100, all sequences with ’N’ will be in your ouput, 
                    set value to 0 if you want no sequences with ”N” in your output""")

parser.add_argument("--n_motifs", required=False, type=int, default=10,
                    help="""Number of motifs able to be founded, DEFAULT = 10""")

parser.add_argument("--size_motifs", required=False, type=int, default=20, help="""Size of the motifs, DEFAULT = 20""")

parser.add_argument("--ncluster", required=False, type=int, default=2,
                    help="Number of clusters able to be calculated, DEFAULT = 2")

args = parser.parse_args()
file_fasta = os.path.abspath(args.input)
alphabet = string.ascii_letters + string.digits
e = datetime.datetime.now()
job_title = '-'.join([e.strftime("%y%m%d_%H-%M-%S"), ''.join(secrets.choice(alphabet).upper() for i in range(10))])
output_job = os.path.join(os.path.abspath(args.output), job_title)
basename_file = 'dataset'


def configure_log():
    try:
        # Loggin configuração coloca o arquivo dentro da pasta log do projeto
        log_file_name = f"{output_job}/{job_title}.log"
        log_format = '[%(asctime)s] [%(name)s]\t%(levelname)s:%(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logging.basicConfig(filename=log_file_name, filemode='w',
                            datefmt=date_format, format=log_format,
                            level=logging.INFO)

        root = logging.getLogger()
        root.setLevel(logging.INFO)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(datefmt=date_format, fmt=log_format))
        root.addHandler(handler)

    except Exception as erro:
        raise erro


def create_tree():
    try:
        path = os.path.join(output_job)
        if not os.path.isdir(path):
            os.mkdir(path)
            os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        else:
            raise "Can't create output folder"

    except Exception as erro:
        raise erro


class Chloe:
    def __init__(self):
        self.main()
        self.len_dataset = int()
        self.len_dataset_after_filter = int()

    def meme_call(self, cluster):
        try:
            module = motif.Motifs_discovery(args.n_motifs, args.size_motifs)
            file = f"{cluster.split('/')[-1]}.fa"
            module.caller(file, cluster)
        except Exception as erro:
            raise erro

    # Função multhreading para o meme
    def multi_meme(self):
        try:
            clusters = [os.path.join(output_job, cluster)
                        for cluster in os.listdir(output_job)
                        if cluster.startswith("cluster_")]

            executor = concurrent.futures.ThreadPoolExecutor(max_workers=int(os.cpu_count()) + 4)
            # perform all tasks in parallel
            futures = [executor.submit(self.meme_call, cluster) for cluster in clusters]
            done, not_done = concurrent.futures.wait(futures, return_when=concurrent.futures.ALL_COMPLETED)

        except Exception as erro:
            raise erro

    # Função controladora da classe do limpador de sequencia
    def cleaner(self):
        try:
            # Inicializando a classe do limpador de sequencia
            module = sequence_cleaner.Sequence_Cleaner(file_fasta,
                                                       output_job,
                                                       args.min_lenght,
                                                       args.max_lenght,
                                                       args.percent_n,
                                                       basename_file)
            # Chamando a funcionalidade cleaner
            module.main()
            self.len_dataset, self.len_dataset_after_filter = module.get_return()
            logging.getLogger("SEQ_CLEANER").info("Total of sequences in input: %s", self.len_dataset)
            logging.getLogger("SEQ_CLEANER").info("Total of sequences after filter: %s", self.len_dataset_after_filter)
        except Exception as erro:
            raise erro

    # Função controladora da classe alinhadora
    def aligner(self, func, name=float('Nan'), file=float("Nan"), output_path=output_job):
        try:
            # Inicializando a classe do alinhador
            module = alignment
            if func == 'main':
                file_input = f"{basename_file}.rmdup.fa"
                # Chamando alinhador de sequencia
                module.main_align(basename_file, f"{output_path}/{file_input}", output_path)

            elif func == 'cluster':
                file_input = os.path.join(output_path, file)
                module.main_align(name, file_input, output_path)

            elif func == 'consensus':
                return module.get_consensus(file)

        except Exception as erro:
            raise erro

    # Função controladora da classe clusterizadora
    def clustering(self):
        try:
            # Inicializando a classe do clusterizador
            module = clustering.Clustering(f"{basename_file}.aligned.fa", output_job, args.ncluster)

            # Chamando a função clusterizadora
            module.main()
        except Exception as erro:
            raise erro

    def blast(self, func, file=float("Nan")):
        try:
            module = blast
            if func == 'createdb':
                module.create_blastdb(file)
        except Exception as erro:
            raise erro

    def get_size_cluster(self,file_path):
        count = 0
        for record in SeqIO.parse(file_path, "fasta"):
            count += 1
        return count
        
    def main(self):
        # Sequence cleaner
        try:
            logging.getLogger("SEQ_CLEANER").info("Calling sequence cleaner")
            self.cleaner()
        except Exception as erro:
            raise erro

        # Alignment
        try:
            logging.getLogger("ALIGNMENT").info("Calling aligner")
            self.aligner('main')
        except Exception as erro:
            raise erro

        # Clustering
        try:
            logging.getLogger("CLUSTERING").info("Calling clusterer")
            self.clustering()
        except Exception as erro:
            raise erro


        # realignment
        try:
            clusters = [dir for dir in os.listdir(output_job) if dir.startswith("cluster_")]
            for cluster in clusters:
                cluster_num = cluster.split("_")[1].split(".")[0]
                file = f"{cluster}.fa"
                name = cluster.split(".")[0]
                logging.getLogger('CLUSTER_ALIGN').info('Align cluster:%s', cluster_num)
                self.aligner('cluster', name, file, os.path.join(output_job, cluster))

        except Exception as erro:
            raise erro

        # creating blastdb
        try:
            clusters = [dir for dir in os.listdir(output_job) if dir.startswith("cluster_")]
            for cluster in clusters:
                file = os.path.join(cluster, f"{cluster}.fa")
                self.blast('createdb', os.path.join(output_job, file))
        except Exception as erro:
            raise erro

        # Multi processing motif discovery
        try:
            logging.getLogger("MEME").info("Calling meme")
            self.multi_meme()
        except Exception as erro:
            raise erro

        # Parsing meme results
        try:
            module = motif
            for cluster in [folder for folder in os.listdir(output_job) if folder.startswith("cluster")]:
                size = self.get_size_cluster(os.path.join(output_job,cluster,f"{cluster}.fa"))
                module.parse_meme(os.path.join(output_job, cluster, 'motif_disc', 'meme.xml'),
                                  os.path.join(output_job, cluster),
                                  size,output_job, cluster)
        except Exception as erro:
            raise erro

        # BLAST
        try:
            clusters = [cluster for cluster in os.listdir(output_job) if cluster.startswith("cluster")]
            for cluster in clusters:
                module = blast.Blast()
                module.search_motif(output_job,
                                    os.path.join(output_job, "MOTIF_QUERY.fa"),
                                    os.path.join(output_job, cluster, f"{cluster}.fa"),
                                    f"{cluster}.blast.txt")
            module.parse_table_blast(output_job)
        except Exception as erro:
            raise erro


if __name__ == "__main__":
    create_tree()
    configure_log()
    logging.getLogger("STARTING").info(f"Starting job...")
    logging.getLogger("STARTING").info(f"Job_ID: {job_title}")
    logging.getLogger("STARTING").info(f"Job_output: {output_job}")
    logging.getLogger("STARTING").info(f"Input fasta: {file_fasta}")
    Chloe()