# Packed Bed Reactor State Estimation

This repository contains the MATLAB scripts and models developed for dynamic modeling and state estimation of a Packed Bed Reactor (PBR) system, as part of my Master's thesis. The project focuses on simulating the reactor dynamics, performing model linearization, implementing PID control, and applying two different Moving Horizon Estimation (MHE) algorithms for state estimation.

# Features

Dynamic modeling of a non-isothermal Packed Bed Reactor (PBR)
Linearization of the dynamic model around operating points
PID controller design and implementation
Moving Horizon Estimation (MHE) with two algorithmic approaches
well-documented MATLAB code

# Tools & Dependencies

MATLAB (R2021b or newer recommended)
CasADi (v3.5.5 or compatible) — for optimization and automatic differentiation

# How to Run the Codes

Instructions:
For the dynamics the main coding file is: catalystDecay_bulk_simulate_parameters.m. Other files were testing different scenarios, proven unsuccesful in thier implementation. For running this file, ensure your path of CasADI is added and you cqn just run that with the helper functions also in your directory.

For the linearization and PID, it is the same story.

# Reference

This repository is part of the work conducted for my Master's thesis at Aalto University. All references are included in the paper; however, the essential references are attached in the repository.
