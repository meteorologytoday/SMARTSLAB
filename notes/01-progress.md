
# TO-DO list

- Test EC-earth and NCAR-CESM model
- Use Kraus-Turner parameterization to see the result


- Understand what processes are actually involved in GFDL ESM2M model
- NCAR-CESM [Download Link](https://www.earthsystemgrid.org/search.html?Project=CMIP5&Experiment=piControl&Ensemble=r1i1p1&Model=CESM1-WACCM&Frequency=Monthly&Product=output1)


# 2019-01-14
- Complete a crude mixed layer model MLMML (Mixed-layer Model - Multiple Layers).
  Convective adjustment is currently only apply on mixed layer. Need to think on applying that to deep ocean.
- For now I do not think it is important to make ocean semi-transparent.

# 2018-11-08
- STAN program applied. Get pretty interesting result. Uncertainty of Q and h are trade-off. Meaning solid physical processes are missing
- Now start to parallelize the STAN program in order to make a map distribution of data.

# 2018-10-24
- Things to finish this week: (1) try using 4096 method. (2) Doing the same for EC-Earth and NCAR-CESM (3) using STAN for a simple case



# 2018-08-29
- Pull out the number to see why places with sea-ice does not do well in simulation. It seems that modifying h so that h = 10m if it gets thinner or negative causes the problem.



# 2018-08-22
- Complete plotting simulation for single points.
- Need to wrap up the notes.

# 2018-08-21
- Transplant to GreenPlanet successfully. Use conda as package manager works fine.
- Complete the simulation for a single point. Need to plot result for different places with RK4/euler scheme.


# 2018-08-20
- Try to transplant the code onto GreenPlanet but still failed. Keep trying.
- Done for the case 1 and 2 and their plottings. For h and Q exist with constant spatial structure, there is unreasonable large uncertainty. For h only with constant spatial strueture, the uncertainty only exists near high latitude coastal regions but h can be negative in a lot of places which is not reasonable.


# 2018-08-16
- Make latex note to include all the methods now i am using.

# 2018-08-14
- Make graph of new equation including dh_dt. Still need to write the latex note down.

# 2018-08-10
- Making plots of 2 methods fitting h and Q. Should clarify these two methods in notes.
- Try to move to green planet. Some package in Julia rely majorly on python. So it is better to setup python first then install julia.

# 2018-08-09
- Did fitting h and Q at the same time. However the result is not good because it seems Q flux is negative generally in boreal summer. This is counter intuitive. What is even weird is that if we assum constant h and fit Q, Q flux will be as expected!


# 2018-08-08
- Write an NetCDFHelper to create nc file more efficiently in Julia.
- Realize the work in 08-07 is actually linear fit. The updated note is [here](https://www.sharelatex.com/read/ffhwmpjxwbht)

# 2018-08-07

- Make the [note](https://www.sharelatex.com/read/ywkvvgyzbmfn) of different approaches.
  

# 2018-08-06

- Make the note of z-average with dynamic boundary


