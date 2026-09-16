#include "check.h"
#include <stdarg.h>

int t_checks = 0;
int t_fails  = 0;
int t_xfail  = 0;
int t_xpass  = 0;

static const char *cur   = "(none)";
static int         cur_f = 0;

void t_begin(const char *name)
{
    cur   = name;
    cur_f = t_fails;
}

void t_end(void)
{
    printf("  %-52s %s\n", cur, (t_fails == cur_f) ? "ok" : "FAIL");
}

void t_fail_at(const char *file, int line, const char *expr, const char *fmt, ...)
{
    t_fails++;
    printf("    FAIL %s:%d\n      %s\n", file, line, expr);
    if (fmt != NULL) {
        va_list ap;
        va_start(ap, fmt);
        printf("      ");
        vprintf(fmt, ap);
        printf("\n");
        va_end(ap);
    }
}

void t_xfail_at(const char *file, int line, const char *id,
                const char *expr, int passed)
{
    t_checks++;
    if (passed) {
        t_xpass++;
        t_fails++;
        printf("    XPASS %s:%d  [%s]\n"
               "      %s\n"
               "      This defect appears FIXED. Remove the EXPECT_FAIL_UNTIL marker.\n",
               file, line, id, expr);
    } else {
        t_xfail++;
        printf("    xfail %s:%d  [%s] (documented defect, still unfixed)\n",
               file, line, id);
    }
}

int t_report(const char *suite)
{
    printf("\n%s: %d checks, %d failures", suite, t_checks, t_fails - t_xpass);
    if (t_xfail) printf(", %d xfail", t_xfail);
    if (t_xpass) printf(", %d XPASS (markers to remove)", t_xpass);
    printf("\n");
    return (t_fails == 0) ? 0 : 1;
}
