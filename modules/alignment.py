#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Created By  : Joao Paulo Ferreira
# Created Date: 15/05/2022
# version ='1.0'
# Module from Chloe: Aligner
# ---------------------------------------------------------------------------
""" Este módulo será responsavel por fazer o alinhamento de sequencia, bem como
funções que estão dentro desse leque."""
# ---------------------------------------------------------------------------
# Imports

import logging
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MafftCommandline
import os

# ------------------------------ALIGNMENT--------------------------------------

class Alignment:

    def get_consensus (self,input_file):
        try:
            logging.getLogger("CONSENSUS").info(f'Opening aligment file: {input_file}')
            alignment = AlignIO.read(open(input_file), 'fasta')
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = summary_align.dumb_consensus(float(0.5))
            return consensus
        except Exception as erro:
            logging.getLogger("CONSENSUS").error(f"Can't generate consensus, erro: {erro}")

    def main(self,name_file,input_file, output_path):
        try:
            logging.getLogger("ALIGNMENT").info("Starting alignment")
            mafft_cline=MafftCommandline(input=input_file, thread = os.cpu_count())
            logging.getLogger("ALIGNMENT").info(f"Sintax: {mafft_cline}")
            stdout, stderr = mafft_cline()
            with open(os.path.join(output_path,f"{name_file}.aligned.fa"), "w") as handle:
                handle.write(stdout)
            logging.getLogger("ALIGNMENT").info("successfully completed!")
        except Exception as erro:
            raise erro