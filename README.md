# Circular elastic RCL example

MATLAB code for the circular truncation experiment in Figure 3. The elastic
field is coupled to phase-extracted P- and S-wave potentials by the weak
interface formulation. Element interiors are eliminated by static
condensation. The circular RCL bilinear form is assembled in radial and
tangential components.

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
only when its parameters and cache version match the current computation.
The data file also records the scaled residual and normwise backward error.

The code was tested with MATLAB R2023a.

## License

This code is distributed under the BSD 3-Clause License. See `LICENSE` for
details.
