# FluidDynTools.jl

This repository contains tools for instruction of fluid dynamics.

[![Build Status](https://github.com/UCLAMAEThreads/FluidDynTools.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/UCLAMAEThreads/FluidDynTools.jl/actions/workflows/CI.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/UCLAMAEThreads/FluidDynTools.jl/HEAD)


To get started with jupyter, please do the following:

1. Download Julia to your computer from https://julialang.org

2. Start a Julia session.

3. Type `]` to enter the Julia package management system. At this prompt, type `add IJulia`. This allows you to run `jupyter`, which is a browser-based computing environment. The 'ju' part stands for Julia; the 'py' part stands for Python. You can return to the Julia prompt by pressing backspace.

If you already have jupyter set up on your computer, start here:

1. Enter the Julia package management system using `]`, if you are not already in it (see above). Type `add https://github.com/UCLAMAEThreads/FluidDynTools.jl`. This will add the FluidDynTools.jl package to your own computer. (In future sessions, you do not need to run this again. However, you should regularly run `update` to make sure you have the most up-to-date version of the package. I.e., type `] update` from the Julia prompt. Please make sure to start a new Julia session after you update, to force a new compilation of the package.)

2. Now press backspace to return to the Julia prompt. At this prompt, type `using FluidDynTools`. This precompiles the package. Note: it might take a few minutes.

3. After you have waited patiently for that to finish, go to a working directory in which you plan to work. You can do this by using `cd("/path/of/your/directory")`.

4. Type `setup_notebooks(pwd())`. This will copy all of the package's Jupyter notebooks into your working directory.

5. Now type `using IJulia` and then `notebook()`. This will start Jupyter in your browser window. (It will download Jupyter via miniconda first if you haven't done it yet.)

6. Navigate to your working directory in the Jupyter file browser. You should see all of your notebooks there. Open `0.0-Index.ipynb`, which will allow you to see the whole set in order.

7. Once you open a notebook, it is recommended that you work on a copy via the menu `File > Make a copy...`. However, if you wish to refresh your notebooks, you can always run the `setup_notebooks` command again, but it will overwrite any of your work.
