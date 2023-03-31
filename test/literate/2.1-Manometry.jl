include("header.jl")

#=
# Manometry
In this notebook, we will learn
* How to use manometers to measure pressure
* How to combine several working fluids in a manometer
* How to calculate for tilted manometers
=#

# ### Set up the module
using MAE103
#-
using Plots

#=
In fluid dynamics, we are often interested in measuring pressure differences
between two points, often because those pressure differences allow us
to determine flow speed.

We'll discuss how we get flow speed from pressure difference later. For now, our
goal is to show how we can measure pressure difference between two
points. The key notion here is that we only need to measure pressure *difference*,
not absolute pressure. Also, one of those "points" might simply be the atmosphere,
where the pressure has some ambient value.

Manometers enable this by turning pressure difference into a measurable height of a
*working fluid* in the manometer. The manometer itself consists of a narrow tube,
bent into a U-shaped form and arranged vertically so that gravity and fluid statics can
be exploited. Sometimes the manometer has additional feature, e.g.,

- one of the channels in the tube might be at an angle rather than vertical
- one or more of the channels might be connected to a "cistern", with a
much wider diameter than the tube.

The working fluid is generally a **liquid**, and we can assume that it does not
change its volume during operation. This fact will be important in the analysis.
=#
#-
#=
### Example 1
We will demonstrate here how we use the fluid statics equation to determine pressure
by reading a measured height. Consider the example setup shown in the figure below:

# <img src="https://raw.githubusercontent.com/UCLAMAEThreads/MAE103/main/notebook/w0022.jpg" alt="manometer" width="250" align="left"/>

We seek to measure the gauge pressure (i.e., relative to ambient) in the air in the tank
by measuring the heights of liquids. Suppose that the oil has specific gravity 0.9. Also,
mercury (Hg) has a very high specific gravity, 13.6: it is a liquid metal.

The heights shown are $h_1 = 36$ inches, $h_2 = 6$ inches, and $h_3 = 9$ inches.
Let's find the gauge pressure in the air, in psi.

First set the values:
=#
h1 = Height(36u"inch")
h2 = Height(6u"inch")
h3 = Height(9u"inch")
γoil = SpecificWeight(SpecificGravity(0.9))
γHg = SpecificWeight(SpecificGravity(13.6))
patm = Pressure(1u"atm");
#=
We work backward from the open end, where pressure is atmospheric. Just to be careful,
we will add this atmospheric pressure, so we are working with absolute pressures,
but then subtract it at the end to find our answer. At point 2, the pressure
increases by $γ_{\mathrm{Hg}}h_3$ from atmospheric:
=#
p2 = Pressure(patm + γHg*h3)
#=
Point 1 on the mercury side is at the same elevation as point 2, and in the same
fluid, so its pressures are equal.
=#
p1 = p2
#=
Also $p_1$ in oil is the same as in mercury. And this pressure is larger than
the pressure in the air by $\gamma_{\mathrm{oil}}(h_1 + h_2)$, so $p_{\mathrm{air}}$
=#
p_air = Pressure(p1 - γoil*(h1+h2))
#=
Remember, this was the absolute pressure. The gauge pressure is obtained by
subtracting the atmospheric pressure:
=#
p_airg = Pressure(p_air-patm)

#=
Finally, subtract the convert to psi
=#
value(p_airg,u"psi")

#=
In summary, the answer was

$$p_{\mathrm{air},g} = γ_{\mathrm{Hg}}h_3 - \gamma_{\mathrm{oil}}(h_1 + h_2)$$
=#

#=
In this example, the manometer used mercury. What if it had used water instead?
What would $h_3$ be (assuming the same gauge pressure in the air)? This is
given by

$$h_3 = \left(p_{\mathrm{air},g} + \gamma_{\mathrm{oil}}(h_1 + h_2)\right)/γ_w$$

where $γ_w$ is the specific weight of water. Let's find out what $h_3$ would be
(in feet).
=#
h3 = Height((p_airg + γoil*(h1+h2))/SpecificWeight(Water))
value(h3,u"ft")

#=
10.2 feet instead of 9 inches! There is big advantage in using a heavy working
liquid in a manometer, because it does not rise as high. It's easier to
get readings.
=#

#=
### Example 2
In the previous example, we pointed out that mercury is useful for translating
a large pressure into a relatively small height change. Sometimes we need the
opposite: if we want to measure a relatively small pressure difference,
we want it to be easy to read on the manometer. For this purpose, we can use an
inclined tube manometer, as in the setup shown in the figure below:

# <img src="https://raw.githubusercontent.com/UCLAMAEThreads/MAE103/main/notebook/w0027.jpg" alt="manometer" width="250" align="left"/>

We read the length $l_2$ along the manometer tube. How does this length come into the
manometer equation? Only the vertical part of $l_2$ enters this equation, i.e.,
$l_2 \sin \theta$

The equation between $A$ and $B$ is

$$p_A + \gamma_1 h_1 - \gamma_2 l_2 \sin \theta - \gamma_3 h_3 = p_B$$

(From $A$, we drop down by $h_1$ with fluid 1, rise by $l_2 \sin \theta$ with fluid 2,
 then rise by $h_3$ with fluid 3.)

Let's assume that fluids 1 and 3 are gases, so their specific weights are negligible.
Then

$$p_B - p_A = \gamma_2 l_2 \sin \theta$$

or, in other words,

$$l_2 = \dfrac{p_B - p_A}{\gamma_2\sin\theta}$$

Let's try this out for a few different angles of inclined tube. We will assume that the
pressure difference is 100 Pa --- a pretty small pressure difference --- and the
working fluid (fluid 2) is water.
=#
dp = Pressure(100)
γ = SpecificWeight(Water)
#=
We'll start with a vertical tube (90 degrees)
=#
Θ = 90π/180
l = Length(dp/γ/sin(Θ))
value(l,u"cm")
#=
Now a tube at 45 degrees
=#
Θ = 45π/180
l = Length(dp/γ/sin(Θ))
value(l,u"cm")
#=
and finally a tube at 10 degrees
=#
Θ = 10π/180
l = Length(dp/γ/sin(Θ))
value(l,u"cm")
#=
Notice how the length goes up as the tube becomes more and more horizontal.
This makes it much easier to read the manometer without being prone to error.
=#
