#!/bin/bash

# Removendo os arquivos .DS_Store recursivamente.
find . -name ".DS_Store" -type f -delete

# Removendo qualquer outro arquivo de metadados do MAC (arquivos que começam com "._") recursivamente.
find . -name "._*" -type f -delete

# Removendo .Trashes and outros que, às vezes, o MAC gera, quando quer.
find . -name ".Trashes" -type d -exec rm -rf {} +

# Remove Python-generated files
find . -type d -name "__pycache__" -exec rm -rf {} +
find . -type d -name "*.egg-info" -exec rm -rf {} +
find . -type d \( -name ".pytest_cache" -o -name ".mypy_cache" -o -name ".ipynb_checkpoints" -o -name ".pytest_cache" \) -exec rm -rf {} +
find . -type f -name "*.py[co]" -delete

# Remove C/C++ compiled files and temporary objects
find . -type f \( -name "*.o" -o -name "*.out" -o -name "*.gch" \) -delete

# (CUIDADO) Removendo todos os arquivos escondidos (i.e. TUDO que começa com ".")
# find . -name ".*" -type f ! -name ".git*" -delete

echo "All tidy mam'!"

