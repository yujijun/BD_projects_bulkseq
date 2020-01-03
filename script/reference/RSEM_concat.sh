#! /bin/bash
# concatenate all result together
paste  *.genes.results | tail -n+2 | cut -f1,5,12,19,26,33,40,47,54,61,68,75,82,89,96,103,110,117,124,131 > ./edgeR.genes.TPM.rsem.txt

paste  *.genes.results | tail -n+2 | cut -f1,6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132 > ./edgeR.genes.TPM.rsem.txt

paste  *.isoforms.results | tail -n+2 | cut -f1,5,13,21,29,37,45,53,61,69,77,85,93,101,109,117,125,133,141,149 > ./edgeR.isoforms.rsem.txt

paste  *.isoforms.results | tail -n+2 | cut -f1,6,14,22,30,38,46,54,62,70,78,86,94,102,110,118,126,134,142,150 > ./edgeR.isoforms.TPM.rsem.txt

