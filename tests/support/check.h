#ifndef PGA_TEST_CHECK_H
#define PGA_TEST_CHECK_H

// Minimal assertion harness.  No dependencies, no framework, no registration
// magic: each test binary has an explicit main() that RUN()s its tests, so the
// control flow reads top to bottom.
//
// MACROS
//   CHECK(expr)            plain assertion
//   CHECK_MSG(expr, ...)   assertion with a printf-style explanation, used
//                          wherever a failure needs context to act on
//   CHECK_EQ_I(a, b)       integer equality; prints got/want on failure
//   RUN(fn)                run one test, print a pass/fail line
//
// EXPECTED FAILURES
//   EXPECT_FAIL_UNTIL(id, expr) marks ONE assertion as a known, documented
//   defect, where id names a catalog entry such as "R-03":
//
//     assertion false -> reported as xfail; the suite stays green
//     assertion true  -> reported as XPASS and the suite FAILS
//
//   The second case is the important one.  When someone fixes the defect the
//   suite goes red with "bug fixed, remove the marker", so a stale marker
//   cannot quietly outlive the bug it documents.  Only the single failing
//   assertion is ever marked; everything else in that test must pass today.
//   A crash, abort or timeout is never an acceptable expected failure.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int t_checks;
extern int t_fails;
extern int t_xfail;
extern int t_xpass;

void t_begin(const char *name);
void t_end(void);
int  t_report(const char *suite);
void t_fail_at(const char *file, int line, const char *expr, const char *fmt, ...);

// Marks a single assertion as a documented, still-unfixed defect.
//   assertion false -> XFAIL, suite stays green
//   assertion true  -> XPASS, suite FAILS ("bug fixed, remove the marker")
void t_xfail_at(const char *file, int line, const char *id,
                const char *expr, int passed);

#define CHECK(expr) \
    do { t_checks++; if (!(expr)) t_fail_at(__FILE__, __LINE__, #expr, NULL); } while (0)

#define CHECK_MSG(expr, ...) \
    do { t_checks++; if (!(expr)) t_fail_at(__FILE__, __LINE__, #expr, __VA_ARGS__); } while (0)

#define CHECK_EQ_I(a, b)                                                      \
    do { long long a_ = (long long)(a), b_ = (long long)(b); t_checks++;      \
         if (a_ != b_)                                                        \
             t_fail_at(__FILE__, __LINE__, #a " == " #b,                      \
                       "got %lld, want %lld", a_, b_);                        \
    } while (0)

#define EXPECT_FAIL_UNTIL(id, expr) \
    t_xfail_at(__FILE__, __LINE__, (id), #expr, (expr) ? 1 : 0)

#define RUN(fn) do { t_begin(#fn); fn(); t_end(); } while (0)

#endif
