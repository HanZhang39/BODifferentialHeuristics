# Biobjective Differential Heuristic
C++ implementations of bi-objective differential heuristics [1].

You can compile the project using the following commands:

``` shell
cmake .
make
```
After typing the above commands. Cmake will generate 4 binaries in the `bin` folder:


+ `bod` runs BOD to compute single-to-all Pareto frontiers.
+ `compress` compresses the Pareto frontier with a given epsilon value
+ `solver_dh` solves bi-objective search problem instances with differential heuristic
+ `solver` is the baseline solver.


You can type `{binary_name} --help` to see the input arguments each binary accepts. Note that the description still needs to be updated.

## Example usage

Please see `example_usage.ipynb` or `run_exp.sh` for example usage

## Reference
[1] Zhang, Han, et al. "Towards effective multi-valued heuristics for bi-objective shortest-path algorithms via differential heuristics." Proceedings of the International Symposium on Combinatorial Search. Vol. 16. No. 1. 2023.
