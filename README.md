# Orbit-Simulation-and-Orbital-transfer
This repository contains multiple functions and scripts to simulate orbit trajectory by solving the Keplerian equation and numerically integrating two body governing equations using custom RK4. In addition, Hohmann transfer and bi-elliptical transfer are included.

# Files

The following MATLAB script and functions are included:

1-  A code that simulates the trajectory of bielliptic transfer between two circular orbits  with given r2: "bi_elliptic.m"
2- A code that converts position and velocity in ECI coordinates to for a given elliptical orbit" ECI2classical.m"
3- A script that computes the trajectory of Hohmann transfer for two circular orbits using RK-4 numerical integration: "Hohmann_trans.m"
4- A script that simulates orbit by numerically integrating the governing equation without/with considering the aerodynamic drag on the spacecraft: "num_orbit_nodrag.m, num_orbit_withdrag.m"
5- A script that simulates a satellite trajectory in cartesian coordinates by solving Kepler's equation" "simulate_orbit.m"
6- A function that solves the Kepler equation to obtain the value of Eccentric anomaly given eccentricity and true anomaly: "solve_Kepler.m"
7- A function that solves the Kepler equation for n orbital periods to find the true anomaly trend: "track_true_anomaly.m"
8-" Main.m,Hohmann.m and Orbit_numerical.m" contain script to solve different orbital mechanics problems using the above functions
