#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "run_cmd.h"

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void sleep_s(double s)
{
    struct timespec ts;
    ts.tv_sec  = (time_t) s;
    ts.tv_nsec = (long) ((s - (double) ts.tv_sec) * 1e9);
    nanosleep(&ts, NULL);
}

// Poll for the child until deadline. Returns 1 and fills *wstatus if reaped,
// 0 if still running at the deadline, -1 if waitpid failed.
static int wait_until(pid_t pid, double deadline, int *wstatus)
{
    while (1)
    {
        pid_t w = waitpid(pid, wstatus, WNOHANG);
        if (w == pid)
            return 1;
        if (w < 0 && errno != EINTR)
            return -1;
        if (now_s() >= deadline)
            return 0;
        sleep_s(0.05);
    }
}

static int open_out(const char *path)
{
    return open(path, O_WRONLY | O_CREAT | O_TRUNC, 0664);
}

RunResult run_with_timeout(char *const argv[], const char *stdout_path,
                           const char *stderr_path, double timeout_s,
                           double kill_grace_s)
{
    RunResult r;
    int out_fd = -1, err_fd = -1, wstatus = 0;
    double t0 = now_s();

    memset(&r, 0, sizeof(r));
    r.status = RUN_SETUP_ERROR;

    if (stdout_path != NULL && (out_fd = open_out(stdout_path)) < 0)
        return r;
    if (stderr_path != NULL && (err_fd = open_out(stderr_path)) < 0)
    {
        if (out_fd >= 0) close(out_fd);
        return r;
    }

    fflush(NULL);
    pid_t pid = fork();
    if (pid < 0)
    {
        if (out_fd >= 0) close(out_fd);
        if (err_fd >= 0) close(err_fd);
        return r;
    }
    if (pid == 0)
    {
        setpgid(0, 0);
        if (out_fd >= 0) { dup2(out_fd, STDOUT_FILENO); close(out_fd); }
        if (err_fd >= 0) { dup2(err_fd, STDERR_FILENO); close(err_fd); }
        execvp(argv[0], argv);
        fprintf(stderr, "run_with_timeout: cannot exec %s: %s\n",
                argv[0], strerror(errno));
        _exit(127);
    }
    setpgid(pid, pid);   // also from the parent, to close the race with kill()
    if (out_fd >= 0) close(out_fd);
    if (err_fd >= 0) close(err_fd);

    double deadline = (timeout_s > 0) ? t0 + timeout_s : 1e300;
    int w = wait_until(pid, deadline, &wstatus);
    if (w < 0)
    {
        r.status = RUN_EXIT_NONZERO;
        r.exit_code = -1;
        r.elapsed_s = now_s() - t0;
        return r;
    }
    if (w == 1)
    {
        r.elapsed_s = now_s() - t0;
        if (WIFEXITED(wstatus))
        {
            r.exit_code = WEXITSTATUS(wstatus);
            r.status = (r.exit_code == 0) ? RUN_OK : RUN_EXIT_NONZERO;
        }
        else if (WIFSIGNALED(wstatus))
        {
            r.signal = WTERMSIG(wstatus);
            r.status = RUN_SIGNALED;
        }
        else
            r.status = RUN_EXIT_NONZERO;
        return r;
    }

    // Over the limit: TERM the group, then KILL, then stop waiting.
    kill(-pid, SIGTERM);
    if (wait_until(pid, now_s() + kill_grace_s, &wstatus) == 0)
    {
        kill(-pid, SIGKILL);
        if (wait_until(pid, now_s() + kill_grace_s, &wstatus) == 0)
        {
            r.status = RUN_STUCK;
            r.elapsed_s = now_s() - t0;
            return r;
        }
    }
    r.status = RUN_TIMEOUT;
    r.elapsed_s = now_s() - t0;
    return r;
}

const char *run_result_describe(const RunResult *r, double timeout_s,
                                char *buf, int buflen)
{
    switch (r->status)
    {
    case RUN_OK:
        snprintf(buf, buflen, "ok (%.1f s)", r->elapsed_s);
        break;
    case RUN_EXIT_NONZERO:
        if (r->exit_code < 0)
            snprintf(buf, buflen, "lost track of the child (waitpid failed)");
        else
            snprintf(buf, buflen, "exited with status %d%s", r->exit_code,
                     r->exit_code == 127 ? " (command not found?)" : "");
        break;
    case RUN_SIGNALED:
        snprintf(buf, buflen, "killed by signal %d", r->signal);
        break;
    case RUN_TIMEOUT:
        snprintf(buf, buflen, "timed out after %.0f s and was killed", timeout_s);
        break;
    case RUN_STUCK:
        snprintf(buf, buflen, "timed out after %.0f s and did not exit after "
                 "SIGKILL (likely stuck in I/O); abandoned", timeout_s);
        break;
    default:
        snprintf(buf, buflen, "could not start (open/fork failed)");
        break;
    }
    return buf;
}

double alntopaf_timeout_from_env(double def_s)
{
    const char *s = getenv("PGA_ALNTOPAF_TIMEOUT");
    char *end;
    double v;

    if (s == NULL || *s == '\0')
        return def_s;
    v = strtod(s, &end);
    if (end == s || v < 0)
        return def_s;
    return v;
}
