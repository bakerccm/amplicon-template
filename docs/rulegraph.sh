#!/bin/bash

# snakemake -j 1 --use-conda --rulegraph | dot -Tpdf >docs/rulegraph.pdf
# snakemake -j 1 --use-conda --rulegraph | dot -Tpng >docs/rulegraph.png

snakemake -j 1 --use-conda --rulegraph >docs/temp
dot -Tpdf docs/temp >docs/rulegraph.pdf
dot -Tpng docs/temp >docs/rulegraph.png
rm docs/temp

# alternatively, use --filegraph or --dag
