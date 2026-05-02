CC     = mpicc
CFLAGS = -O3 -Wall -Wno-unused-result -fno-strict-aliasing -I./lib -I./src
LIBS   = -lpthread -lm -lz

LIB_SRCS = lib/GDB.c lib/gene_core.c lib/libfastk.c lib/ONElib.c \
           lib/MSDsort.c lib/align.c lib/alncode.c lib/ANO.c

CORE_SRCS   = src/extract.c src/pack.c src/tables.c src/work_split.c src/record.c
PHASE1_SRCS = src/pga-extract.c $(CORE_SRCS)
PHASE2_SRCS = src/pga-sort.c    $(CORE_SRCS) src/sort.c

all: distributedpanga pga-extract pga-sort

distributedpanga: src/main.c src/ingest.c $(LIB_SRCS)
	$(CC) $(CFLAGS) -o distributedpanga src/main.c src/ingest.c $(LIB_SRCS) $(LIBS)

# Phase 1 single-process extractor. Uses gcc rather than mpicc since it has no MPI dependency.
pga-extract: $(PHASE1_SRCS) $(LIB_SRCS)
	gcc $(CFLAGS) -o pga-extract $(PHASE1_SRCS) $(LIB_SRCS) $(LIBS)

# Phase 2 single-process extractor + msd_sort.
pga-sort: $(PHASE2_SRCS) $(LIB_SRCS)
	gcc $(CFLAGS) -o pga-sort $(PHASE2_SRCS) $(LIB_SRCS) $(LIBS)

verify-phase1: pga-extract
	bash scripts/verify_phase1.sh

verify-phase2: pga-sort
	bash scripts/verify_phase2.sh

clean:
	rm -f distributedpanga pga-extract pga-sort
