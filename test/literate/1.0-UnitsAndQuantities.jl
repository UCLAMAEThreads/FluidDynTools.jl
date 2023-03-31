include("header.jl")

#=
# Fluid properties, flow quantities, and units
In this notebook, we will discuss *fluid properties* and *flow quantities*,
and the systems of units we use for them. First, some basic definitions:
- Fluid properties are properties of the material (the fluid)
- Flow quantities are characteristics of the flow of this fluid
We will also introduce some syntax we will use in this notebook and those that follow for dealing with units.
=#

# ### Set up the module
using MAE103

#=
We will generally focus on SI units, and these will be the default system
for displaying quantities. However, we also need to be familiar
with Imperial (sometimes called "English") units, which arise in many
situations in engineering. In the examples below, we will show that the tools
in these notebooks allow us to easily convert from one system to another.

For any quantity, we can see what the default units are by using `default_unit`.
For example, for pressure,
=#
default_unit(Pressure)
#-
#=
### Fluid properties
Let's start by discussing the basic properties of a fluid.

#### Density
The **density** provides a measure of the amount of fluid material per unit volume.
It is measured in units of mass/volume:
=#
default_unit(Density)

#=
For example, the density of water at a reference temperature of 15.6 degrees C is
=#
Density(Water)

#=
and for air at temperature of 15 degrees C and pressure of 1 atmosphere is
=#
Density(Air)

#=
Notice that water is around 800 times denser than air. This fact is very important
in fluid mechanics!
=#

#=
Remember that the density of gases, like air, depend on pressure and temperature (by
the ideal gas law), so the density of a gas may be sensitive to the local conditions.
This will happen in flows traveling relatively fast ([Q: compared to what?]). In such a case, density is a flow quantity,
not a fluid property, and we have to determine its value as part of the problem.

However, when a gas is traveling relatively slow [Q: compared to what?], we can
often treat a gas as having constant density, so it can be treated as a fluid property.

And liquid density is generally not sensitive. In fact, usually we can assume that the
density in a liquid is constant and uniform. We can treat
density in a liquid as a fluid property rather than a flow quantity.
=#
#-
#=
Also note that the density of seawater is larger (at average salt concentration)
than that of freshwater
=#
Density(Seawater)


#=
What if we want to set the density in other units? For example, in imperial
units, we would usually set it with lbm/ft^3. (In the notebooks, we use `lbm` for
pound (mass), and `lbf` is used for pound (force).) To set the units of a quantity,
rather than rely on the default units, we use the syntax to follow the number with
 `u"units"`. It will automatically convert it to the default units. For example, 2 lb/ft^3:
=#
Density(2u"lbm/ft^3")

#=
and to report the value of a quantity in other units, we use the `value`
function. For example, to get the density of water in lb/ft^3:
=#
value(Density(Water),u"lbm/ft^3")



#=
#### Viscosity
The internal friction in the fluid, called the *viscosity* (or, more specifically, the *dynamic viscosity*), given
by the Greek symbol $\mu$. [We can get this symbol by typing `\mu+TAB`.] The
viscosity controls the relationship between *shear stress*  $\tau$ (the frictional
force per unit area) and the *strain rate*, given by the gradient of velocity,
$\mathrm{d}u/\mathrm{d}y$. This latter quantity measures the difference in speeds in adjacent layers
of fluid, and larger differences suggest more shear stress. Viscosity is
the proportionality constant

$$\tau = \mu \dfrac{\mathrm{d}u}{\mathrm{d}y}$$

We will learn
much more about viscosity later, but for now, it is sufficient to
know that viscosity has units of kg/m/s:
=#
μw = Viscosity(Water)
#-
μa = Viscosity(Air)

#=
Note that water is much more viscous the air:
=#
μw/μa
#=
This is probably intuitive to you.
=#
#=
Note that both of these viscosities are much smaller than that of glycerin:
=#
Viscosity(Glycerin)
#-
#=
We will also occasionally make use of the ratio between viscosity and density.
This is called the *kinematic viscosity*. We use the symbol $\nu$ for this
[obtained by typing `\mu+TAB`.]
=#
νw = KinematicViscosity(Water)
#-
νa = KinematicViscosity(Air)



#=
#### Surface tension
The surface tension is a property associated with liquid interfaces, and particularly,
liquid interfaces with gases. The surface tension has units of force per unit length,
because if we imagine "cutting" a bit of the interface from the rest of it,
then this cut would form a perimeter of the snipped part of the interface, and
surface tension would act along this perimeter, representing how much the rest
of the liquid interface was pulling on it.

Surface tension also happens to have units of energy per unit area, so it
is sometimes referred to as "surface energy". Different liquids have different
surface energies, depending on the strength of their inter-molecular forces.
=#
SurfaceTension(Water)
#-
SurfaceTension(Glycerin)

#=
### Flow quantities
Now, let's discuss quantities that describe the fluid flow, or
at least, the *state* of the fluid. It is important to understand that
these quantities are, in general, *field quantities*: they vary from location
to location, and perhaps vary over time. So each of them should be
thought of as a *function* of the spatial coordinates, $(x,y,z)$ and time $t$.
Finding these functions is often our ultimate goal in solving a problem.
=#

#=
#### Pressure
Pressure represents the average force that the molecules exert per unit
area of surface. It is important to understand that pressure acts the same in every
direction. That is, it is an *isotropic* quantity. This means that, no matter
what the orientation of the surface, the pressure acts the same on it, and further,
it only acts *perpendicular* to the surface.

The default SI unit of pressure is the *Pascal* (Pa, equal to 1 N/m^2, or 1 kg/m/s^2).
=#
Pressure(20)

#=
But there are many other units for pressure in use. Some common ones. The
atmosphere (atm) represents the standard ambient pressure.
=#
p = Pressure(1u"atm")
#=
So 1 atm is about 101325 Pa, or 101.325 kPa. Let's see this written in other units
=#
value(p,u"psi")
#-
value(p,u"mmHg")
#=
These are all *absolute* pressures. For many of the flows we will
study, the absolute pressure will not be important. Consider,
for example, the flow of water through a pipe. Only the *difference* between
the pressure at the inlet and the pressure at the outlet matters for
driving this flow; the absolute pressure is irrelevant. So
we can also define a `PressureDifference`, measured in the
same units:
=#
PressureDifference(50)

#=
Pressure is a flow quantity. In general, we need to find its values
as part of the problem.
=#

#=
#### Velocity
The velocity of the fluid describes how fast it is moving in each direction,
and has units of length/time.
=#
u = Velocity(20)
#-
#=
In other units,
=#
value(u,u"cm/s")
#-
value(u,u"ft/s")

#=
We know that velocity is a vector, so it has, in general, three components $(u,v,w)$.
It is also a field quantity. This means that each of the three components
depends on $(x,y,z)$ and $t$. That's a lot of detail to determine, but once
we determine it, we know *everything* about a fluid flow.
=#
