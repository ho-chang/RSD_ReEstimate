# Parameter Search for Hybrid Algebraic Attacks on the RSD problem

This README provides instructions on how to run and use the implementation files for searching the parameters $(\eta_{1}, \eta, \eta_{\rm pre})$ related to the proposed algorithms presented in the paper:

**“Improving Hybrid Algebraic Attacks on the Regular Syndrome Decoding Problem with Balanced Guessing Strategy”**

## File Description

This repository contains the following files:

- `findPar.c`
  - `findPar.c` is a C source file used to compile a parameter-search program.
  - The compiled executable supports the search for optimal complexity estimates across various instances of the Regular Syndrome Decoding (RSD) problem, in accordance with the algorithms and parameter choices discussed in the paper.
  - `findPar_Beaver.c` is a C source file used to compile a parameter-search program for the PCGs Beaver [22].

- `exp.py`
  - `exp.py` is a Python script provided to facilitate the reproduction of the experimental results presented in the paper, including those reported in Table 1 .
  - It enables users to execute the parameter configurations considered in the paper and to compare the resulting complexity estimates with the reported values.

- `paper_attack.py`
  - `paper_attack.py` is a Python script that provides a small-scale implementation of the proposed attack, intended to illustrate its overall procedure and computational behavior on instances with ($n \leq 20$ and $\beta = n/h \in \mathbb{Z}$).
.


## Usage & Example

To find the parameter tuple $(\eta_{1}, \eta, \eta_{\rm pre})$ for an instance of $\textup{RSD}(n,k,h)$, first compile `findPar.c` as follows:
```
gcc -O3 -march=native -flto -fopenmp -fno-math-errno -std=c11 findPar.c -lm -o rsd
```
This command generates the executable rsd.

For example, to search for the minimum complexity for the instance $\textup{RSD}(1024,652,106)$ over $\mathbb{F}_{2^{128}}$, as reported in Table 1 with $\xi(N)=2$, run:
```bash
./rsd 1024 652 106 0 2
```
or
```bash
OMP_NUM_THREADS=16 ./rsd 1024 652 106 0 2
```
It then outputs the resulting complexity, the corresponding parameter tuple $(\eta_{1}, \eta, \eta_{\rm pre})$, and other relevant information.
Note that, for the binary field, an example with $a=1.0$ (while the default value is set to $0.5$) is given as follows:

```bash
./rsd 1024 652 106 1 1.0
```
It searchs for the minimum complexity for the instance $\textup{RSD}(1024,652,106)$ over $\mathbb{F}_{2}$ with $\xi(N)=1.0 \cdot N$.


---

The script `exp.py` can be used to reproduce the experimental results reported in the paper.
To display the available options, run:
```bash
python ./exp.py
```
The script provides 14 options for reproducing the data presented in the paper.

For example, to reproduce the results for the parameter group indexed by 0, run:
```bash
python ./exp.py 0
```
This command searches for the parameters corresponding to group 0.

--- 
The following commands provide an example of how to execute paper_attack.py with a small parameter set.

Windows PowerShell:

```powershell
python paper_attack.py 18 10 3 --eta1 2 --eta2 2 --pre 2 --seed 2 --outer-budget 20 --inner-budget 20 --trace
```

Linux:

```bash
python3 paper_attack.py 18 10 3 --eta1 2 --eta2 2 --pre 2 --seed 2 --outer-budget 20 --inner-budget 20 --trace
```

In this example, the positional arguments 18 10 3 specify the problem parameters, while `--eta1`, `--eta2`, and `--pre` determine the corresponding algorithmic parameters.
The `--seed` option fixes the random seed for reproducibility, and `--outer-budget` and `--inner-budget` specify the search budgets used in the outer and inner stages, respectively.
The `--trace` option enables detailed output of the intermediate steps of the algorithm.


## Notes on Reducing the Search Time (for optimal parameter $(\eta_{1}, \eta, \eta_{\rm pre})$)

We provide two additional remarks concerning the running time of the parameter search.

First, when the given RSD instance involves very large parameters, for example when ($n > 2^{20}$), an exhaustive search for optimal parameters may require a substantial amount of time.
This issue is particularly relevant for the PCG Beaver instances, where a full exhaustive search can be computationally expensive.
However, for the PCG Beaver parameters considered in this work, the corresponding quadratic systems are already solvable by linearization.
Therefore, we provide the file `findPar_Beaver.c`, which searches for the minimum complexity only under the first zero guessing.

Second, in practice, one may obtain good parameters for a given RSD instance within a significantly shorter running time by restricting the search space.
We observed that many optimal parameter choices satisfy the condition $n_{2} = n_{\rm pre} + 1$.
Although we do not provide a proof that this condition always returns an optimal parameter choice, imposing this restriction substantially reduces the number of candidates considered in the exhaustive search.
Users who wish to perform a faster, restricted search may modify `findPar.c` accordingly.
