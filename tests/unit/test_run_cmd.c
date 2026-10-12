// ============================================================================
// R-05 -- External converter run with a time limit.  See tests/CATALOG.md.
// ============================================================================
//
// THE PROBLEM THIS GUARDS.  pga-mpi-merge and pga-merge used to run
// `ALNtoPAF ... 2>/dev/null` through system(): no time limit, stderr thrown
// away, and a failure printed but ignored (exit 0).  In job 59686567 the
// converter hung for over two hours after the pipeline had finished, and
// there was no log to say why.  run_with_timeout() (src/run_cmd.c) replaces
// system(); these tests drive it with small shell commands instead of
// ALNtoPAF so every outcome can be produced on purpose.
//
// Not covered: RUN_STUCK (a child that survives SIGKILL). That needs a process
// in uninterruptible I/O, which a unit test cannot create on demand.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#include "check.h"
#include "run_cmd.h"

static char tmpdir[64];

static void path_in_tmp(char *out, size_t n, const char *name)
{
    snprintf(out, n, "%s/%s", tmpdir, name);
}

static int file_equals(const char *path, const char *want)
{
    char buf[256];
    FILE *f = fopen(path, "r");
    size_t got;
    if (f == NULL) return 0;
    got = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[got] = '\0';
    return strcmp(buf, want) == 0;
}

static RunResult run_sh(const char *script, double limit, double grace)
{
    char out[128], err[128];
    char *argv[] = { "sh", "-c", (char *) script, NULL };
    path_in_tmp(out, sizeof(out), "out");
    path_in_tmp(err, sizeof(err), "err");
    return run_with_timeout(argv, out, err, limit, grace);
}

static void r05_ok_captures_streams(void)
{
    char out[128], err[128];
    RunResult r = run_sh("echo hello; echo oops >&2", 10, 1);
    path_in_tmp(out, sizeof(out), "out");
    path_in_tmp(err, sizeof(err), "err");
    CHECK_EQ_I(r.status, RUN_OK);
    CHECK_EQ_I(r.exit_code, 0);
    CHECK_MSG(file_equals(out, "hello\n"), "stdout not captured");
    CHECK_MSG(file_equals(err, "oops\n"), "stderr not captured");
}

static void r05_nonzero_exit(void)
{
    RunResult r = run_sh("exit 3", 10, 1);
    CHECK_EQ_I(r.status, RUN_EXIT_NONZERO);
    CHECK_EQ_I(r.exit_code, 3);
}

static void r05_command_not_found(void)
{
    char out[128], err[128];
    char *argv[] = { "pga-no-such-command-r05", NULL };
    path_in_tmp(out, sizeof(out), "out");
    path_in_tmp(err, sizeof(err), "err");
    RunResult r = run_with_timeout(argv, out, err, 10, 1);
    CHECK_EQ_I(r.status, RUN_EXIT_NONZERO);
    CHECK_EQ_I(r.exit_code, 127);
}

static void r05_killed_by_signal(void)
{
    RunResult r = run_sh("kill -SEGV $$", 10, 1);
    CHECK_EQ_I(r.status, RUN_SIGNALED);
    CHECK_EQ_I(r.signal, SIGSEGV);
}

static void r05_timeout_terminates(void)
{
    RunResult r = run_sh("sleep 30", 0.5, 2);
    CHECK_EQ_I(r.status, RUN_TIMEOUT);
    CHECK_MSG(r.elapsed_s < 5, "took %.1f s to stop a 0.5 s limit", r.elapsed_s);
}

// Ignores SIGTERM, so only the SIGKILL escalation can stop it.
static void r05_timeout_escalates_to_kill(void)
{
    RunResult r = run_sh("trap '' TERM; sleep 30", 0.5, 0.5);
    CHECK_EQ_I(r.status, RUN_TIMEOUT);
    CHECK_MSG(r.elapsed_s < 5, "took %.1f s", r.elapsed_s);
}

// The limit applies to the whole process group, including grandchildren.
static void r05_timeout_kills_group(void)
{
    RunResult r = run_sh("sleep 30 & wait", 0.5, 2);
    CHECK_EQ_I(r.status, RUN_TIMEOUT);
    CHECK_MSG(r.elapsed_s < 5, "took %.1f s", r.elapsed_s);
}

static void r05_unwritable_output(void)
{
    char *argv[] = { "true", NULL };
    RunResult r = run_with_timeout(argv, "/nonexistent-dir-r05/out", NULL, 10, 1);
    CHECK_EQ_I(r.status, RUN_SETUP_ERROR);
}

static void r05_timeout_env(void)
{
    unsetenv("PGA_ALNTOPAF_TIMEOUT");
    CHECK(alntopaf_timeout_from_env(600) == 600);
    setenv("PGA_ALNTOPAF_TIMEOUT", "45", 1);
    CHECK(alntopaf_timeout_from_env(600) == 45);
    setenv("PGA_ALNTOPAF_TIMEOUT", "0", 1);
    CHECK(alntopaf_timeout_from_env(600) == 0);      // 0 = no limit
    setenv("PGA_ALNTOPAF_TIMEOUT", "abc", 1);
    CHECK(alntopaf_timeout_from_env(600) == 600);
    setenv("PGA_ALNTOPAF_TIMEOUT", "-5", 1);
    CHECK(alntopaf_timeout_from_env(600) == 600);
    unsetenv("PGA_ALNTOPAF_TIMEOUT");
}

int main(void)
{
    char cmd[128];
    snprintf(tmpdir, sizeof(tmpdir), "/tmp/test_run_cmd.XXXXXX");
    if (mkdtemp(tmpdir) == NULL)
    {
        perror("mkdtemp");
        return 1;
    }
    printf("test_run_cmd\n");
    RUN(r05_ok_captures_streams);
    RUN(r05_nonzero_exit);
    RUN(r05_command_not_found);
    RUN(r05_killed_by_signal);
    RUN(r05_timeout_terminates);
    RUN(r05_timeout_escalates_to_kill);
    RUN(r05_timeout_kills_group);
    RUN(r05_unwritable_output);
    RUN(r05_timeout_env);
    snprintf(cmd, sizeof(cmd), "rm -rf %s", tmpdir);
    if (system(cmd) != 0)
        fprintf(stderr, "warning: could not remove %s\n", tmpdir);
    return t_report("test_run_cmd");
}
