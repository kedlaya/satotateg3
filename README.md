# satotateg3
This repository consists of code associated to the paper [*Sato-Tate groups of abelian threefolds*](http://arxiv.org/abs/2106.13759
) (Fite, Kedlaya, Sutherland).
This code depends on the following software packages:

- SageMath (tested with version 9.3)
- Magma (tested with version 2.25-5)

Unless otherwise noted, the relevant package can be inferred from the file extension (`.sage` for SageMath, `.m` for Magma).

The code is organized into the following directories.

- Section5:
    - `make_groups.sage`: construct presentations of the groups of type N.
    - `compute_group_inclusions.sage`: compute inclusions of groups of type N (Proposition 5.39).
    - `beukers_smyth.sage`: run the Beukers-Smyth algorithm (Remark 5.40).
    - `check_abelian_group.sage`: check Table 4, Proposition 5.23, Proposition 5.25 (Remark 5.42).
    - `compute_group_ids.sage`: compute isomorphism classes of the component groups of type N (Remark 5.43).
- Section6:
    - See separate README file in this folder.
- Section7:
    - `GSp6Largeimage.m`: test whether the Jacobian of a given curve of genus 2 or 3 has large Galois image based on its Euler factors at small primes.
- Section8:
    - `check_galois_groups.sage`: check various assertions about explicit curves and their Jacobians (Examples 8.4--8.34). This code requires Magma to be installed.
    - `check_twisting_representations.sage`: check groups presentations used for twisting powers of CM elliptic curves (Examples 8.17--8.34).
    - `Example*.m`: compute L-functions for Examples 8.17--8.34.
    - `common.m`: some utility functions used in `Example*.m`.
    - `*FrobeniusContext.txt`: some precomputed data used in `Example*.m`.  These files will automatically be recomputed if absent, but this takes much longer.
    - The `moments` directory contains the file `mzdata.txt` which contains moment statistics and point densities for each Sato-Tate group (as computed by the code in the Section6 folder), and the files `moments_*.txt` contain numerically computed moment statistics and point densities for the base changes of Examples 8.4--8.34 to every subfield of their endomorphism fields.  The bash script `checkmoments.sh` compares this data. and reports any variance greater than one percent.
