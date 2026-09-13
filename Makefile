CC     = mpicc
CFLAGS = -O3 -Wall -Wno-unused-result -fno-strict-aliasing -I./lib -I./src
LIBS   = -lpthread -lm -lz

# Build with byte-resolution LCPs by default.
# This matches FastGA's stock GIXmake (which is also compiled with -DLCPs);
# stripping it produces a 0/1 boundary marker that is NOT what FastGA's
# downstream merge expects.
#
# Both lib/MSDsort.c (msd_sort writes real LCP values) and src/dsort.c
# (cross-rank LCP boundary fixup uses compute_kmer_lcp) see the same flag.
#
# Opt out (rarely useful) with:  make LCPS=0 pga-mpi
LCPS ?= 1
ifeq ($(LCPS),1)
CFLAGS += -DLCPs
endif

LIB_SRCS = lib/GDB.c lib/gene_core.c lib/libfastk.c lib/ONElib.c \
           lib/MSDsort.c lib/align.c lib/alncode.c lib/ANO.c lib/RSDsort.c

CORE_SRCS   = src/extract.c src/pack.c src/tables.c src/work_split.c src/record.c
PHASE1_SRCS = src/pga-extract.c $(CORE_SRCS)
PHASE2_SRCS = src/pga-sort.c    $(CORE_SRCS) src/sort.c
PHASE3_STAGE1_SRCS = src/pga-mpi-stage1.c $(CORE_SRCS) src/sort.c src/dsort.c
PHASE3_SRCS        = src/pga-mpi.c        $(CORE_SRCS) src/sort.c src/dsort.c

# Phase 4 single-rank prototype: extract + sort + adapter + FastGA merge/align/PAF.
# No MPI — does not pull in dsort.c.
PHASE4_SRCS        = src/pga-merge.c $(CORE_SRCS) src/sort.c \
                     src/kmer_adapter.c src/fastga_pipeline.c

all: pga-extract pga-sort pga-mpi-stage1 pga-mpi pga-merge pga-mpi-merge

# Phase 1 single-process extractor. Uses gcc rather than mpicc since it has no MPI dependency.
pga-extract: $(PHASE1_SRCS) $(LIB_SRCS)
	gcc $(CFLAGS) -o pga-extract $(PHASE1_SRCS) $(LIB_SRCS) $(LIBS)

# Phase 2 single-process extractor + msd_sort.
pga-sort: $(PHASE2_SRCS) $(LIB_SRCS)
	gcc $(CFLAGS) -o pga-sort $(PHASE2_SRCS) $(LIB_SRCS) $(LIBS)

# Phase 3 stage 1 — gather-to-rank-0 baseline. Throwaway; deleted after stage 2 lands.
pga-mpi-stage1: $(PHASE3_STAGE1_SRCS) $(LIB_SRCS)
	$(CC) $(CFLAGS) -o pga-mpi-stage1 $(PHASE3_STAGE1_SRCS) $(LIB_SRCS) $(LIBS)

# Phase 3 stage 2 — real distributed pipeline (extract + Alltoallv + per-rank sort).
pga-mpi: $(PHASE3_SRCS) $(LIB_SRCS)
	$(CC) $(CFLAGS) -o pga-mpi $(PHASE3_SRCS) $(LIB_SRCS) $(LIBS)

# Phase 4 — single-rank adapter prototype. gcc, no MPI.
pga-merge: $(PHASE4_SRCS) $(LIB_SRCS)
	gcc $(CFLAGS) -o pga-merge $(PHASE4_SRCS) $(LIB_SRCS) $(LIBS)

# Phase 4b — distributed merge + alignment. mpicc, links pga-mpi's stage 2
# (dsort.c) and pga-merge's adapter + lifted FastGA pipeline.
PHASE4B_SRCS = src/pga-mpi-merge.c $(CORE_SRCS) src/sort.c src/dsort.c \
               src/kmer_adapter.c src/fastga_pipeline.c src/contig_assignment.c \
               src/seed_exchange.c
pga-mpi-merge: $(PHASE4B_SRCS) $(LIB_SRCS)
	$(CC) $(CFLAGS) -o pga-mpi-merge $(PHASE4B_SRCS) $(LIB_SRCS) $(LIBS)

verify-phase1: pga-extract
	bash scripts/verify_phase1.sh

verify-phase2: pga-sort
	bash scripts/verify_phase2.sh

verify-phase3-stage1: pga-mpi-stage1
	bash scripts/verify_phase3_stage1.sh

verify-phase3-stage2: pga-mpi
	bash scripts/verify_phase3_stage2.sh

# Phase 3 stage 2 LCP regression. Builds with LCPS=1 and runs Check A
# (within-group LCP=KMER), Check B (between-group LCP matches local
# computation), plus a cross-check against GIXshow.
verify-phase3-stage2-lcp:
	bash scripts/verify_phase3_stage2_lcp.sh

clean:
	rm -f pga-extract pga-sort pga-mpi-stage1 pga-mpi pga-merge pga-mpi-merge

# ---------------------------------------------------------------------------
# Tests.  See tests/CATALOG.md for the catalog (test IDs are cited in the sources).
#
# Tests compile with the production CFLAGS so -DLCPs matches the LCPS=1
# default; a test built without it silently exercises the 0/1 LCP branch.
# Each binary links only the production sources it actually needs.
# ---------------------------------------------------------------------------
TEST_CFLAGS  = $(CFLAGS) -I./tests/support
TEST_SUPPORT = tests/support/check.c tests/support/synth_gdb.c \
               tests/support/ref_syncmer.c tests/support/rec_build.c
TEST_GDBLIB  = lib/GDB.c lib/gene_core.c lib/ONElib.c lib/ANO.c
TEST_EXSRCS  = src/extract.c src/pack.c src/tables.c src/work_split.c src/record.c

test-extract:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -o tests/bin/test_extract tests/unit/test_extract.c \
	    $(TEST_SUPPORT) $(TEST_EXSRCS) $(TEST_GDBLIB) $(LIBS)
	tests/bin/test_extract

# Same sources, reduced scan window: exercises the multi-window re-scan path
# on a 200-base contig instead of a 20 Mbp one.
test-extract-window:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -DSCAN_MAX=64 -o tests/bin/test_extract_window \
	    tests/unit/test_extract_window.c \
	    $(TEST_SUPPORT) $(TEST_EXSRCS) $(TEST_GDBLIB) $(LIBS)
	tests/bin/test_extract_window

test-lcp-sort:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -o tests/bin/test_lcp_sort tests/unit/test_lcp_sort.c \
	    $(TEST_SUPPORT) src/sort.c $(TEST_EXSRCS) $(TEST_GDBLIB) lib/MSDsort.c $(LIBS)
	tests/bin/test_lcp_sort

test-record-split:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -o tests/bin/test_record_split tests/unit/test_record_split.c \
	    tests/support/check.c src/record.c src/work_split.c $(LIBS)
	tests/bin/test_record_split

test: test-record-split test-extract test-extract-window test-lcp-sort \
      test-dsort test-assign test-adapter test-exchange-plan

# Slow paths (production allocation growth) plus the sanitizers.
test-asan:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -fsanitize=address,undefined -O1 -g \
	    -o tests/bin/test_lcp_sort_asan tests/unit/test_lcp_sort.c \
	    $(TEST_SUPPORT) src/sort.c $(TEST_EXSRCS) $(TEST_GDBLIB) lib/MSDsort.c $(LIBS)
	tests/bin/test_lcp_sort_asan --slow
	gcc $(TEST_CFLAGS) -fsanitize=address,undefined -O1 -g \
	    -o tests/bin/test_extract_asan tests/unit/test_extract.c \
	    $(TEST_SUPPORT) $(TEST_EXSRCS) $(TEST_GDBLIB) $(LIBS)
	tests/bin/test_extract_asan

test-clean:
	rm -rf tests/bin

# dsort.h includes mpi.h for its prototypes, so this builds with $(CC); the
# functions under test are pure and the binary runs as a single process.
test-dsort:
	@mkdir -p tests/bin
	$(CC) $(TEST_CFLAGS) -o tests/bin/test_dsort tests/unit/test_dsort.c \
	    $(TEST_SUPPORT) src/dsort.c src/sort.c $(TEST_EXSRCS) $(TEST_GDBLIB) \
	    lib/MSDsort.c $(LIBS)
	tests/bin/test_dsort

test-assign:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -o tests/bin/test_assign tests/unit/test_assign.c \
	    tests/support/check.c src/contig_assignment.c $(LIBS)
	tests/bin/test_assign

test-adapter:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -o tests/bin/test_adapter tests/unit/test_adapter.c \
	    tests/support/check.c tests/support/rec_build.c tests/support/ref_syncmer.c \
	    tests/support/synth_gdb.c src/kmer_adapter.c src/record.c src/tables.c \
	    $(TEST_GDBLIB) $(LIBS)
	tests/bin/test_adapter

test-exchange-plan:
	@mkdir -p tests/bin
	gcc $(TEST_CFLAGS) -o tests/bin/test_exchange_plan tests/unit/test_exchange_plan.c \
	    tests/support/check.c src/seed_exchange.c $(LIBS)
	tests/bin/test_exchange_plan

test-mpi:
	@mkdir -p tests/bin
	$(CC) $(TEST_CFLAGS) -o tests/bin/test_stage2 tests/mpi/test_stage2.c \
	    tests/support/check.c tests/support/rec_build.c tests/support/ref_syncmer.c \
	    tests/support/synth_gdb.c src/dsort.c src/sort.c $(TEST_EXSRCS) \
	    $(TEST_GDBLIB) lib/MSDsort.c $(LIBS)
	bash tests/run_mpi.sh
