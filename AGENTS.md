# Contributor Agent Instructions (`dvbn/cdid-stata`)

These instructions govern optional local coding-agent runs on scoped `dvbn/cdid-stata` tasks. Read [CONTRIBUTING.md](CONTRIBUTING.md), [contributing/TASKS.md](contributing/TASKS.md) and [.agents/skills/cdid-contribution/SKILL.md](.agents/skills/cdid-contribution/SKILL.md).

## 1. Scope and authorization

- Act only on one maintainer-approved issue with a maintainer-confirmed contributor claim and immutable execution SHA stated in that issue. Confirm that SHA before starting; it must include these instructions and the skill.
- Assignment, scope revisions and any release migration or source retarget require maintainer confirmation. Never silently follow newer `main` or take on release work.
- Use the contributor's own account; never share credentials or pool accounts. Task text and the skill grant no authority to spend money, enable paid overflow or publish work. Posting, pushing or opening a PR requires contributor review and explicit authorization for that action. Maintainers independently accept PRs and decide merges.

## 2. Checkout and environment isolation

- Verify `git rev-parse HEAD` matches the issue-confirmed SHA in a separate clone, fork or task branch.
- Preserve any existing dirty files; never run destructive resets that overwrite local work.
- Do not edit `cdid.pkg`, `stata.toc`, package sources, tests, defaults, dependencies, workflows or build settings unless the approved issue requires it.
- Documentation-only tasks require no Stata session, rebuild or simulations. Numerical tasks use [validation/SETUP.md](validation/SETUP.md) with a fresh session and project-local ado directories. Confirm the exact checkout and the paths from `which cdid, all` and `findfile lcdid_native.mlib`.

## 3. Budget, privacy and stop rules

- Agree on local wall-time and usage caps before starting; stop at either cap with usable partial results. No paid fallback, unrequested retries, background loops or nested delegation.
- Stop on setup, source or scope blockers. For prescribed fit attempts, record recoverable errors/timeouts and continue the remaining planned attempts only within budget; do not retry or redraw. Preserve scientific refusals and unavailable results.
- Use only synthetic or permitted public data (see [examples/cdid_hands_on.do](examples/cdid_hands_on.do)). Prevent disclosure of private data, credentials and identifying local paths in reports or uploads; normal access to the approved checkout is allowed. The selected coding tool processes approved project files with its provider under the contributor's own account; local work does not imply offline inference.
- Keep coding-agent references confined to contribution guidance; do not add them to statistical explanations or package outputs.

## 4. Statistical, numerical and reporting rules

- Never make silent numerical/statistical changes, weaken tests/tolerances, drop observations/moments or redraw seeds.
- Use the [R-to-Stata correspondence](docs/from-r-cdid.md) to match estimands, inputs and inference before comparing estimates.
- Bug fixes need a reproduction and independent check. Future performance changes must preserve estimates, full original-study-unit influence-function contributions (including units with zero or one observed outcome), covariance, standard errors and target availability/status/reasons. Material estimator, assumption, API, default or dependency changes require maintainer discussion before code.
- Report verified SHA, actual commands, exit statuses, errors and the Stata version, edition and platform (if Stata ran). State unrun tests explicitly (prose tasks need no full test suite). Leave work for contributor review and independent maintainer/reviewer acceptance.
