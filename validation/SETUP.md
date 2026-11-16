# Checking a contribution

Record the commit with `git rev-parse HEAD`. From a clean Stata session, use the [project installation instructions](../docs/installation.md#install-into-a-project-folder) with fresh ado directories. For a local checkout, use its absolute path in `net install cdid, from(...)` and the checkout's `csdid` subdirectory for the optional ordinary-mode dependency. Check `which cdid, all`, `findfile lcdid_native.mlib` and, for ordinary mode, `csdid version`. If your task changes Mata source, rebuild in a fresh Stata session using [the source build instructions](../docs/installation.md) before testing.

Run a focused synthetic example and the relevant portions of the [hands-on do-file](../examples/cdid_hands_on.do). Record the Stata version, edition, operating system, commands, return codes and any checks you could not run. Compare results against an independent fixture when changing estimates or inference. The licensed runtime is required for these checks.

Public GitHub Actions verifies the distribution manifest, install catalogs and local Markdown links. Run the same check locally with `python3 validation/check_distribution.py`. For a change to a distributed file, update its manifest entry using the final byte count and SHA-256. For compiled changes, update the native library entry as well after rebuilding and running the numerical checks.
