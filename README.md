# chloe
C.H.L.O.E. - Conserve High-throughput Localization Of genomic Elements
---
CHLOE is a tool for clustering and searching conserved regions in viral genomes. Its purpose is to perform a cross-search of conserved regions among the formed clusters to reduce errors caused by low sequence representativeness in the dataset.

## Key Features

- **Advanced Clustering**: CHLOE employs k-medoids clustering algorithms to group related viral sequences, making it easier to identify genetic patterns and variations.

- **Conserved Region Search**: The tool conducts a thorough search within each cluster to identify highly conserved regions, which can play a critical role in viral evolution studies and vaccine development. For this, the MEME tool is used.

## Why CHLOE?

- **Accelerate Your Research**: With the ability to efficiently group sequences and identify conserved regions, CHLOE saves time and resources, speeding up progress in viral research.

- **Enhance Accuracy**: Eliminate the uncertainty caused by low-representative sequences, obtaining more accurate and reliable results in your genomic analyses.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)

## Requirements

### Hardware Requirements
- Threads = 8
- RAM = 8GB for small datasets (Maximum 2000 sequences)
- Storage = 128GB SSD

### Software Requirements
- Python v3.10.12
- Docker v23.0.1
- MAFFT v7.490

### Docker images
- memesuite/memesuite v5.5.0
- ncbi/blast v2.13.0

## Installation

```bash
- apt-get update
- apt-get upgrade
- git clone https://github.com/biojpferreira/chloe
- docker pull memesuite/memesuite:5.5.0
- docker pull ncbi/blast:2.13.0
- pip install -r requirements.txt
```

---

## Getting Started

- **To use the CHLOE virtual environment**: source environment/bin/activate
- **To view all possible configurations**: python chloe.py -h

## Options
| Argument | full_name | Description |
|----------|-----------|-------------|
| -h | --help | show this help message and exit |

---

  --input INPUT, -i INPUT
                        Path for fasta file.
  --output OUTPUT, -o OUTPUT
                        Path job output.
  --min_lenght MIN_LENGHT
                        the user defines the minimum length. DEFAULT = 0, it means you don’t have to care about the minimum length
  --max_lenght MAX_LENGHT
                        the user defines the maximum length. DEFAULT = 0, it means you don’t have to care about the maximum length
  --percent_n PERCENT_N
                        The user defines the percent of N is allowed. DEFAULT = 100, all sequences with ’N’ will be in your ouput, set value
                        to 0 if you want no sequences with ”N” in your output
  --n_motifs N_MOTIFS   Number of motifs able to be founded, DEFAULT = 10
  --size_motifs SIZE_MOTIFS
                        Size of the motifs, DEFAULT = 20
  --ncluster NCLUSTER   Number of clusters able to be calculated, DEFAULT = 2
