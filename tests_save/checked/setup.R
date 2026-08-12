# Runs once, before any test file.
#
# The tetrachoric paths (measure = "r"/"z" from 2x2 tables, phi, or chi-square)
# delegate to metafor::escalc(measure = "RTET"/"ZTET"), which in turn requires
# the 'mvtnorm' package. metafor only loads mvtnorm lazily, the first time a
# tetrachoric correlation is actually computed.
#
# In a full-suite run many packages are pulled in via `pkg::` calls before the
# tetrachoric tests execute. If mvtnorm's namespace has not been loaded by then,
# `requireNamespace("mvtnorm")` can intermittently fail (e.g. once R's
# "maximal number of DLLs" ceiling is under pressure). metafor then aborts with
# "Please install the 'mvtnorm' package to compute this measure.", and
# metaConvert's own tetrachoric helper (.tet_r) silently returns NA.
#
# Loading mvtnorm here, at the very start of the session while few DLLs are in
# use, makes its namespace stay resident for the whole run so every later
# requireNamespace("mvtnorm") succeeds. This removes the order-dependent
# tetrachoric failures without altering any effect-size logic.
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  loadNamespace("mvtnorm")
}
