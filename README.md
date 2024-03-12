# Description

This is a solver written to solve the two-dimensional incompressible Navier-Stokes equations discretized using the Finite-Volume Method on a staggered grid [[1]](#1). The problem solved is lid-driven cavity flow with the upper wall provided a velocity and the other walls kept stationary.

![Lid Driven Cavity](/pictures/problem.png)

## Running the code

- Install the Fortran 90 compiler. 
- Run the bash file `run.sh`. It compiles the code using the `lblas` and `llapack` flags, runs the executable, and deletes it.

```sh
sh run.sh
```

## Results

The solution is post-processed using MATLAB. The streamlines for different values of Reynold's number are as follows:

![Streamlines](/pictures/stream.png)

The code was validated against ANSYS Fluent. The profiles of horizontal and vertical velocities are plotted along the centre lines.

![Validation](/pictures/val.png)

The following figure shows the comparison of Gauss-Seidel and TDMA solvers:

![Performance](/pictures/perf.png)



## References
<a id="1">[1]</a> 
Versteeg, Henk Kaarle, and Weeratunge Malalasekera. An introduction to computational fluid dynamics: the finite volume method. Pearson education, 2007.
