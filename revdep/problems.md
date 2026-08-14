*Wow, no problems at all. :)*

<!--
metaConvert 2.0.1, checked 2026-08-14.

metaumbrella 1.1.0 was checked against this build with the development metaConvert
installed into an isolated library (see revdep/check-revdeps.R). Result: Status 1 NOTE,
which is metaumbrella's OWN pre-existing unused-Imports note for 'pwr'. Its examples,
both test files and its vignettes all pass.

The previous contents of this file reported

    metaumbrella -- Newly broken
    checking package dependencies ... ERROR
    Package required but not available: 'metaConvert'

That was NOT a regression. It was revdepcheck failing to install metaConvert into its
own private library -- the same Windows failure mode ("package is in use and will not be
installed", and file.rename ... Access denied) that makes revdepcheck unusable here. A
missing dependency in the checking environment is indistinguishable in the log from a
package that genuinely broke, which is exactly how that entry came to be recorded
against metaConvert 1.0.3 and then carried forward unread.

Use revdep/check-revdeps.R instead; it verifies each reverse dependency's own hard
dependencies are resolvable BEFORE checking, and reports a local gap as a local gap.
-->
