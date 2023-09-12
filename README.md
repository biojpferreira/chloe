# chloe
C.H.L.O.E. - Conserve High-throughput Localization Of genomic Elements
---
CHLOE é uma ferramenta de custerização e busca de regiões conservadas em genomas virais. Seu propósito é fazer uma busca cruzada de regiões conservadas entre os clusteres formados para diminuir os erros de baixa representatividade de sequencias no dataset.
Uma breve descrição do seu projeto. Inclua o propósito, a funcionalidade e talvez até mesmo uma captura de tela ou um logotipo.

## Principais Recursos

- **Clusterização Avançada**: O CHLOE emprega algoritmos de clusterização kmedoids para agrupar sequências virais relacionadas, tornando mais fácil identificar padrões e variações genéticas.

- **Busca de Regiões Conservadas**: A ferramenta realiza uma busca minuciosa em cada cluster para identificar regiões altamente conservadas, que podem desempenhar um papel crítico em estudos de evolução viral e desenvolvimento de vacinas para isso é utilizado a ferramenta MEME.

## Por que CHLOE?

- **Acelere Suas Pesquisas**: Com a capacidade de agrupar sequências e identificar regiões conservadas de maneira eficiente, o CHLOE economiza tempo e recursos, acelerando o progresso em pesquisas virais.

- **Aumente a Precisão**: Elimine a incerteza causada por sequências de baixa representatividade, obtendo resultados mais precisos e confiáveis em suas análises genômicas.

## Índice

- [Requisitos](#requisitos)
- [Instalação](#instalação)
- [Como Usar](#como-usar)
- [Contribuição](#contribuição)
- [Licença](#licença)

## Requisitos

### Requisitos de hardware
- Threads = 8
- RAM = 8Gb para datasets pequenos (Maximo 2000 sequências)
- Armazenamento = 128Gb SSD

### Requisitos de software
- Docker v23.0.1
- mafft v7.490

### Docker images
- memesuite/memesuite v5.5.0
- ncbi/blast v2.13.0

## Instalação

```bash
apt-get update
apt-get upgrade

git clone https://github.com/biojpferreira/chloe
docker pull memesuite/memesuite:5.5.0
docker pull ncbi/blast:2.13.0
pip install -r requirements.txt
---
## Primeiros passos
