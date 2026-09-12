// E-11: multi-window extraction.
//
// Built with -DSCAN_MAX=64 (see PR-1 in tests/CATALOG.md), so a 200-base contig
// forces at least three scan windows and exercises the re-window pointer
// arithmetic at src/extract.c:196-203, `Get_Contig_Piece(...) - beg`.
//
// The result is compared against the independent oracle -- not merely against
// the single-window run -- so a seam bug that corrupted both paths identically
// could not pass. test_extract runs these same fixtures at the production
// SCAN_MAX, so agreement between the two scan paths follows transitively.

#include <stdio.h>
#include "check.h"
#include "extract_case.h"
#include "tables.h"

static void e11_multi_window(void)
{
    int i;
    CHECK_MSG(SCAN_MAX < WINDOW_FIXTURE_LEN,
              "fixture must exceed SCAN_MAX (%d) to force a re-window", SCAN_MAX);
    printf("      [SCAN_MAX=%d, fixture len=%d -> >= %d windows]\n",
           SCAN_MAX, WINDOW_FIXTURE_LEN, (WINDOW_FIXTURE_LEN + SCAN_MAX - 1) / SCAN_MAX);
    for (i = 0; i < WINDOW_FIXTURE_N; i++) {
        char label[64];
        snprintf(label, sizeof label, "E-11 SCAN_MAX=%d seed=%u",
                 SCAN_MAX, WINDOW_FIXTURE_SEEDS[i]);
        extract_case_vs_oracle(WINDOW_FIXTURE_LEN, WINDOW_FIXTURE_SEEDS[i], label);
    }
}

// A longer contig crosses many more seams.
static void e11_many_windows(void)
{
    extract_case_vs_oracle(1000, 8675309u, "E-11 long, many windows");
}

int main(void)
{
    init_comp_table();
    printf("test_extract_window (SCAN_MAX=%d)\n", SCAN_MAX);
    RUN(e11_multi_window);
    RUN(e11_many_windows);
    return t_report("test_extract_window");
}
