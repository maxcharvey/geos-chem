# Startup history traversal regression

Run `bash test/implementation/run_history_startup_test.sh` with GNU Fortran.
The test compiles the complete `History_Write` routine extracted from the
current source and supplies minimal state/I/O interfaces. It does not compile
or run GEOS-Chem, load model inputs, or require NetCDF.

The regression covers all 24 traversal permutations of BoundaryConditions,
SolarCache, Processes and Restart. It checks one mapped startup boundary
write, no other startup writes, unchanged chemistry state, subsequent due
writes and alarm advancement, restart file closing, absent/empty/BC-only
lists, not-due writes, and invalid startup units.

For a negative control, set `HISTORY_WRITE_SOURCE` to an unpatched
`History/history_mod.F90`. The original routine fails when BoundaryConditions
is not the first visited collection.

I/O and `History_Update` are stubbed. Passing these tests proves traversal and
the exercised writer branches, not actual NetCDF timestamps, diagnostic
accumulation neutrality, precision, or model restart equality. Those require
the separate one-hour compiled-model qualification.
