Relativistic thin disk model
==============================

Relativistic thin disk model by Novikov & Thorne (1973).

This model implements the formulae presented in 
Novikov&Thorne(1973) and Page&Thorne(1974) papers.

* https://ui.adsabs.harvard.edu/abs/1974ApJ...191..499P
* https://ui.adsabs.harvard.edu/abs/1973blho.conf..343N


### Compilation

```bash
$ make
```

After compilation, two files are produced: 
* a shared object library `disk-nt.so`
* a test executable `disk-nt-test`

The shared object library can be used in custom programs. The test 
executable demonstrates a simple dump of quantities for a given 
choice of spin and mass accretion rate. See the `main()` function in `disk-nt.c`.

