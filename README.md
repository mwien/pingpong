# pingpong

This is a solver for the one-sided crossing minimization problem (OCM) as a submission to the PACE 2024 Exact Track. The input format and I/O behavior is as specified by the [PACE](https://pacechallenge.org/2024/).

The tool was developed by *Max Bannach, Florian Chudigiewitsch, Kim-Manuel Klein* and *Marcel Wienöbst*.

# Algorithm

The solver first reduces the problem to a directed weighted feedback arc set (FAS) instance as, e.g., described in [1]. This instance is then solved by an Implicit Hitting Set formulation as given in [2,3]. We refine this approach in multiple ways. Most notably, we make use of multiple heuristics to add hitting set constraints: heuristics for the hitting set problem to find a good solution for the current formulation and heuristics for FAS to find cycles, which are not hit and thus have to be added as constraints. Inspired by the interplay of these heuristics, we name our solver ``pingpong''.

Correctness is guaranteed as the program will only terminate if an ILP solver finds an optimal solution for the hitting set instance that is checked to eliminate all cycles in the FAS instance. If this is not the case, then further constraints are added and the procedure repeats. 

Disclaimer: The run-time per instance ca vary between runs even on the same hardware as the procedure of adding constraints uses randomness. 

1. Alexander Dobler: *[A Note on the Complexity of One-Sided Crossing Minimization of Trees](https://arxiv.org/abs/2306.15339).* (Technical Report, 2023)
2. Martin Grötschel, Michael Jünger, and Gerhard Reinelt: A cutting plane algorithm for the linear ordering problem. Operations Research 32 (1984).
3. Ali Baharev, Hermann Schichl, Arnold Neumaier, and Tobias Achterberg: An exact method for the minimum feedback arc set problem. ACM Journal of Experimental Algorithmics 26 (2021).

# Dependencies
The following open source [crates](https://crates.io) are used. They are automatically downloaded and compiled when the solver is build using *Cargo*. 
- [bitset](https://crates.io/crates/bit-set)
- [highs](https://crates.io/crates/highs)
- [itertools](https://crates.io/crates/itertools)
- [rand](https://crates.io/crates/rand)

# Build
pingpong is implemented in [Rust](https://www.rust-lang.org) and can simply be build using [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html):

```
cargo build --release
```

# Run
After the build is completed, the tool can either be executed directly via

```
./target/release/pingpong < <instance.gr>
```

or by using Cargo

```
cargo run --release < <instance.gr>
```

