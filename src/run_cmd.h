#ifndef DPGA_RUN_CMD_H
#define DPGA_RUN_CMD_H

// Run an external command with a wall-clock limit (replaces system() for the
// ALNtoPAF conversion). No shell is involved: argv[0] is looked up on PATH.

typedef enum {
    RUN_OK = 0,        // exited with status 0
    RUN_EXIT_NONZERO,  // exited with a non-zero status (127: could not exec)
    RUN_SIGNALED,      // killed by a signal it did not get from us
    RUN_TIMEOUT,       // exceeded the limit; we sent SIGTERM then SIGKILL
    RUN_STUCK,         // still alive after SIGKILL (e.g. uninterruptible I/O);
                       // abandoned, so the caller never hangs on it
    RUN_SETUP_ERROR    // could not open the output files or fork
} RunStatus;

typedef struct {
    RunStatus status;
    int       exit_code;  // valid for RUN_OK / RUN_EXIT_NONZERO
    int       signal;     // valid for RUN_SIGNALED
    double    elapsed_s;
} RunResult;

// stdout_path / stderr_path are truncated and receive the child's output
// (NULL leaves the stream inherited). timeout_s <= 0 means no limit.
// The child runs in its own process group; on timeout the whole group gets
// SIGTERM, then SIGKILL after kill_grace_s, then we wait up to kill_grace_s
// more before giving up with RUN_STUCK.
RunResult run_with_timeout(char *const argv[], const char *stdout_path,
                           const char *stderr_path, double timeout_s,
                           double kill_grace_s);

// Short human-readable description, e.g. "timed out after 600 s".
// Returns a pointer into buf.
const char *run_result_describe(const RunResult *r, double timeout_s,
                                char *buf, int buflen);

// PGA_ALNTOPAF_TIMEOUT (seconds) if set and parseable, else def_s.
double alntopaf_timeout_from_env(double def_s);

#endif
