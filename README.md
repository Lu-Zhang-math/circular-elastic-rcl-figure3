# Circular elastic RCL example

MATLAB code for the circular truncation experiment in Figure 3. The elastic
field is coupled to phase-extracted P- and S-wave potentials by the weak
interface formulation. The modal RCL blocks use the positive test weight
`gamma = r/tau`, giving a Petrov--Galerkin discretization with trial space
`Q_h` and test space `gamma Q_h`. On the circular interface `gamma = 1`, and
on the outer boundary `gamma m_b beta_b = 1`.

Run in MATLAB:

```matlab
results = run_figure3;
```

The computation uses `d = [0.01, 0.1, 1]`, `omega = [5, 10, 50]`, and
`L = 2:14`. There are ten angular elements and five radial layers: four in
the physical domain and one in the RCL. The polynomial degrees
`(N1,N2,N3)` are `(14,40,24)`, `(18,45,30)`, and `(30,60,40)` for
`omega = 5`, `10`, and `50`, respectively. They remain fixed for every
value of `d` and `L`.

Results are written to `results`. Completed parameter cases are stored in
`results/cases`, so an interrupted run can be resumed. A cached case is used
only when its parameters, cache version, and source hash match the current
computation. Invalid or incomplete cache files are recomputed. The result
files record the errors, residuals, fitted rates, fit diagnostics, MATLAB
version, and source hash.

Validate a completed run with:

```matlab
report = validate_figure3;
```

The code was tested with MATLAB R2023a.

## License

This code is distributed under the BSD 3-Clause License. See `LICENSE` for
details.
