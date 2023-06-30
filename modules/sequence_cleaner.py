#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Created By  : Joao Paulo Ferreira
# Created Date: 15/05/2022
# version ='1.0'
# Module from Chloe: Sequence Cleaner
# ---------------------------------------------------------------------------
""" Este módulo será responsavel por limpar sequencia, com parametro de:
comprimento, duplicatas, porcentagem de caracteres ambiguos"""
# ---------------------------------------------------------------------------
# Imports

import logging
import re
from Bio import SeqIO

# -----------------------------SEQUENCE CLEANER-------------------------------


class Sequence_Cleaner:

    def __init__(self, fasta_file, output, min_lenght, max_lenght, percent_n, basename_file):
        self.output = output  # Caminho de output forecido pelo usuário
        self.fasta_file = fasta_file  # Arquivo fasta contendo as sequencias a serem analisadas
        self.max_lenght = max_lenght  # Tamanho máximo de sequencia permitido
        self.min_lenght = min_lenght  # Tamanho minimo de sequencia permitido
        self.percent_n = percent_n  # Porcentagem máximo de caracteres ambiguous permitidos
        self.total_input_seqs = int()  # numero de sequencias fornecidas pelo usuario
        self.total_used_seqs = int()  # numero de sequencia utilizáveis pós filtro
        self.basename_file = basename_file  # Nome do dataset

    def main(self):
        logging.getLogger("SEQ_CLEANER").info('Initializing the Sequence Set Cleaner')
        self.clean_sequence()
        self.remove_dup()
        self.translate_dataset()

    def get_return(self):
        return self.total_input_seqs, self.total_used_seqs

    # TODO: Mudar a forma de tradução da sequencia
    def translate_seq(self, seq):
        if float.is_integer(len(str(seq)) / 3):
            return seq.translate()

        elif float.is_integer((len(str(seq)) - 1) / 3):
            return seq[:-1].translate()

        else:
            return seq[:-2].translate()

    def clean_sequence(self):
        try:
            # Cria dicionário para receber as sequencias
            sequences = {}
            pattern = re.compile(r"[^ACGT]")
            # Using the Biopython fasta parse we can read our fasta input
            for seq_record in SeqIO.parse(self.fasta_file, "fasta"):

                # Take the current sequence
                sequence = str(seq_record.seq).upper()
                self.total_input_seqs = self.total_input_seqs + 1

                # Rename ID
                seq_record.id = str(seq_record.id).replace(' ', '-')
                # Check if the current sequence is according to the user parameters
                if (
                        int(self.min_lenght) <= len(sequence) <= int(self.max_lenght)
                        and (float(sequence.count("N")) / float(len(sequence))) * 100 <= self.percent_n
                        and (float(sequence.count("-")) / float(len(sequence))) * 100 <= self.percent_n
                ):
                    self.total_used_seqs = self.total_used_seqs + 1
                    sequence = pattern.sub("", sequence)
                    # If the sequence passed in the test "is it clean?" and it isn't in the
                    # hash table, the sequence and its id are going to be in the hash
                    if sequence not in sequences:
                        sequences[sequence] = seq_record.id
                    # If it is already in the hash table, we're just concatenate the ID
                    # of the current sequence to another one that is already in the hash table
                    else:
                        sequences[sequence] += "_" + seq_record.id

            # Write the clean sequences

            # Create a file in the same directory where you ran this script
            with open(f"{self.output}/{self.basename_file}.fa", "w+") as output_file:
                # Just read the hash table and write on the file as a fasta format
                for sequence in sequences:
                    output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
            logging.getLogger("SEQ_CLEANER").info("The sequence has been cleaned!")

        except Exception as erro:
            logging.getLogger("SEQ_CLEANER").error("Can't clean the input file.")
            raise erro

    def remove_dup(self):
        try:
            seen = set()
            records = []

            for record in SeqIO.parse(f"{self.output}/{self.basename_file}.fa", "fasta"):
                if record.seq not in seen:
                    seen.add(record.seq)
                    records.append(record)

            SeqIO.write(records, f"{self.output}/{self.basename_file}.rmdup.fa", "fasta")
            logging.getLogger("SEQ_CLEANER").info("Duplicates has been removed!")

        except Exception as erro:
            logging.getLogger("SEQ_CLEANER").error("Can't remove duplicates from input file.")
            raise erro

    def translate_dataset(self):
        try:
            records = []

            for record in SeqIO.parse(f"{self.output}/{self.basename_file}.fa", "fasta"):
                record.seq = self.translate_seq(record.seq)
                records.append(record)

            SeqIO.write(records, f"{self.output}/{self.basename_file}.rmdup.fa", "fasta")
            logging.getLogger("TRANSLATE").info("Sequences has been translated!")

        except Exception as erro:
            logging.getLogger("TRANSLATE").error("Can't translate sequences from input file.")
            raise erro
