# Notes on different derivations
- [Derivation of exact MLD energy](https://www.sharelatex.com/read/xqqjhbzxnqvb)
- [Mixed-layer Depth and( Q-flux Numerical Calculation](https://www.sharelatex.com/read/ffhwmpjxwbht)


# Notes on grid transformation
- NCL : [rcm2rgrid] (https://www.ncl.ucar.edu/Document/Functions/Built-in/rcm2rgrid.shtml)


# Notes on different Models

## GFDL-ESM2G
 - Mixed-Layer calculation: `parameterizations/vertical/GOLD_mixed_layer.F90`. It uses the Kraus-Turner-like bulk mixed layer based on the work [Simulation of the Atlantic Circulation with a Coupled Sea Ice-Mixed Layer-Isopycnal General Circulation Model. Part I: Model Description](https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1993)023%3C0808:SOTACW%3E2.0.CO;2)
 - Grid points: search `NXTOT` and `NYTOT`.




# Resolustion
- ESM2G vertical resolution: 59 density layers beneath the 2 mixed layers and 2 buffer layers. [ESM2G technical paper](https://journals.ametsoc.org/doi/pdf/10.1175/JCLI-D-11-00560.1), [Southern Ocean Circulation and Eddy Compensation in CMIP5 Models](https://journals.ametsoc.org/doi/pdf/10.1175/JCLI-D-12-00504.1)



# Notes on work

## Linear Regression

- Constant h and non-constant h in shallow water linear regression may result in completely different inversion. Please run the code `07_analytic_test/plot.jl`. Through this test we know that wrong assumption to the model is critical.

## Newtonian Fitting Formulation and Result ( [Link](https://www.overleaf.com/read/nhtyxqvdgjvr) )
