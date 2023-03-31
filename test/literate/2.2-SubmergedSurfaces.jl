include("header.jl")

#=
# Submerged surfaces
In this notebook, we will learn how to calculate force and moment on
flat surfaces immersed in static fluids
=#

# ### Set up the module
using MAE103
#-
using Plots


#=
# Example
A circular gate of 3 m diameter is located in an inclined wall (60 degrees
from horizontal) of a large reservoir filled with water. A horizontal shaft runs through
the center of the gate. The shaft is at a depth of 12 m. The other side
of the gate is exposed to air at atmospheric pressure.
* What is the force on the gate?
* Where does this force effectively act?
* What is the moment applied by the water on the gate shaft?
=#
#-
#=
Here, we are working with a simple, symmetric shape (a circle), and we
know that the centroid of a circle is at its center. Since the shaft runs
through the center, we know that the centroid is at a depth of 12 m.
=#
hc = Depth(12)
d = Diameter(3)
Θ = 60*π/180

#=
The force is easy to compute, since it is just the specific weight times
area times depth.

$$F_R = \gamma h_c A$$

(The atmospheric pressure contributes nothing in this problem, since it
acts the same on both sides.)
=#
A = Area(Circle,d/2)
FR = Force(SpecificWeight(Water)*A*hc)

#=
To find the center of pressure, we use the formula

$$y_R = \dfrac{I_{xc}}{y_c A} + y_c$$

This the location along the $y$ axis, which runs from the free surface downward, parallel
to the plate. In the formula, $I_{xc}$ is the second area moment about the $x$ axis running
through the centroid of the plate, and $y_c$ is the $y$ coordinate of the centroid.

But $y_c$ is simply $h_c/\sin\theta$:
=#
yc = Length(hc/sin(Θ))

#=
and $I_{xc}$ we can compute with a handy function:
=#
Ixc = SecondAreaMomentX(Circle,d/2)

#=
Now put them together to get the center of pressure
=#
yR = Length(Ixc/(yc*A) + yc)

#=
How far is this from the centroid?
=#
yR - yc

#=
Remember, it is lower than the centroid, because there are higher pressures
acting below the centroid than above it.

To determine the moment $M_R$ about the shaft, we simply need to remember that the
force is acting at the center of pressure, so the moment is equal to the force times the distance
of this center of pressure from the shaft (i.e., from the centroid)
=#
MR = Moment((yR - yc)*FR)

#=
#### A follow-up question
What happens to these values as the plate gets deeper?

Obviously force gets larger. What about $y_R$ and $M_R$? Let's set the depth
at 20 m, and recalculate:
=#
hc = Depth(20)
FR = Force(SpecificWeight(Water)*A*hc)
yc = Length(hc/sin(Θ))
yR = Length(Ixc/(yc*A) + yc)
yR - yc
#=
so it gets a little closer to the centroid as the plate gets lower. What about
the moment?
=#
MR = Moment((yR - yc)*FR)
#=
It's exactly the same! In fact, **the moment about the centroid does not
depend on how deep the plate is**.
=#

#=
#### Another follow-up question
What happens if the plate becomes more vertical? We'll set $h_c$ to 12 m again
and make $Θ = 90^\circ$. Only $y_c$ changes:
=#
hc = Depth(12)
Θ = 90*π/180
FR = Force(SpecificWeight(Water)*A*hc)
#=
Note that **force is not affected** by a change in angle, since only depth
of the centroid matters.
=#
yc = Length(hc/sin(Θ))
yR = Length(Ixc/(yc*A) + yc)
yR - yc

#=
So the center of pressure is a little further now from the centroid than
when it was at 60 degrees.

The moment:
=#
MR = Moment((yR - yc)*FR)

#=
So **moment gets larger as the surface becomes more vertical**. This is because the
center of pressure has moved away from the centroid. This is intuitive if you
think about the other extreme: a completely horizontal surface. In that situation,
the pressure acts uniformly across the plate, so there is **zero net moment**
on a horizontal surface.
=#
