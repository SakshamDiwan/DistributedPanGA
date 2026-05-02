#!/bin/bash
cd ~/DistributedPanGA/notes
gcc -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing -DLCPs \
    -I ~/FASTGA \
    -o GIXmake_dump GIXmake_dump.c \
    ~/FASTGA/MSDsort.c ~/FASTGA/libfastk.c ~/FASTGA/ONElib.c \
    ~/FASTGA/ANO.c ~/FASTGA/GDB.c ~/FASTGA/gene_core.c \
    -lpthread -lm -lz
echo "Build done: $?"
