---
layout: default
title: Home
nav_exclude: true
---

# Quill: 3D QED particle-in-cell code

**Quill** (simulator for **QU**antum effects in **I**ntense **L**aser-p**L**asma  interactions) is a fully three-dimensional parallel particle-in-cell (PIC) code developed at the [Institute of Applied Physics RAS](https://ipfran.ru/), Nizhny Novgorod, Russia.
To our knowledge, it was the first PIC code with implementation of the Monte Carlo QED approach to investigate the development of electron–positron cascades.

The code is able to model the following processes using the Monte Carlo technique:

* photon emission by an electron in the strong field, with radiation reaction effects;
* electron–positron pair creation from gamma photons (Breit–Wheeler process).

The Maxwell solvers implemented in the code are FDTD, NDFX (the scheme used in A. Pukhov's VLPL code), and hybrid five-point FDTD (the scheme reduces numerical Cherenkov instability).
The particles pushers implemented in the code use Vay or Boris scheme.

To learn more about Quill, read the full documentation on:
* [Dependencies and how to build Quill](build)
* [Input files and how to run Quill](run)
* [Analyzing the results](analyze)