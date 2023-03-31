include("header.jl")

#=
# Example calculation in a pipe system
In this notebook we will see an example of how to solve a problem in a pipe
system, using the 1-d steady-flow energy equation.

In our calculations, we will use the units that we discussed in notebook 1.0
in order to take advantage of the useful tools that come with them.
=#

# ### Set up the module
using MAE103

#=
## Example: Find the flow rate
In this example (Example 8.10 in the textbook), we are to determine the change of
flow rate that occurs after we modify a design for a fume hood. The fume hood
involves a fan that pumps air and noxious gas out of an enclosed region,
through a duct, and out to the atmosphere where it can safely mix with the
ambient air. The existing design involves a fan with a short duct of diameter
8 inches, and the losses correpond to a loss coefficient of $K_L = 5$. This
leads to a flow rate of $\dot{Q}_0 = 9$ ft$^3$/s, which is safely within
the regulated range of 6 to 12 ft$^3$/s.

The new design involves a long 100 ft pipe of galvanized iron (leading to
major loss from viscosity) and a new total loss coefficient of $K_L = 10$.
We wish to see if the new design's volume flow rate falls within the
regulated range. We will assume that the fan's head does not change in
the redesign.
=#
#-
#=
First, let's put in the parameters. We treat the gas as air, since any noxious
gas species are assumed to be at low concentrations.

Galvanized iron has a typical roughness of 0.0005 ft.
=#
Q0 = VolumeFlowRate(9u"ft^3/s")
KL0 = 5
D = Diameter(8u"inch")
μ = Viscosity(Air)
ρ = Density(Air)
ϵ = Height(0.0005u"ft")
g = Gravity()
#=
Now let's calculate the flow velocity in the original design
=#
A = Area(π/4*D^2)
V0 = Velocity(Q0/A)
value(V0,u"ft/s")

#=
And now, we can calculate the fan head $h_p$, since this is equal
to the kinetic energy of the exiting flow plus the head loss. (Pressures
are equal and ambient in the fume hood and at the exit, and elevation
changes are negligible.

$$h_p = \dfrac{V^2}{2g} + K_L\dfrac{V^2}{2g}$$
=#
hp = Head((KL0+1)*V0^2/(2*g))
value(hp,u"ft")

#=
Now we will analyze the modified system. We add 100 ft of duct after the
fan and increase the loss coefficient to 10:
=#
L = Length(100u"ft")
KL = 10.0

#=
The energy equation now becomes

$$h_p = \dfrac{V^2}{2g} + K_L\dfrac{V^2}{2g} + \dfrac{fL}{D}\dfrac{V^2}{2g}$$

The unknown in this equation is $V$, but $f$ (the friction factor) is also
unknown. Let's solve for $V$, pretending that we know $f$:

$$V = \left(\dfrac{2gh_p}{1 + K_L + fL/D}\right)^{1/2}$$

To get $f$, we need Reynolds number and roughness coefficient.
And to get Reynolds number, we need $V$. So, we clearly must iterate.

Our approach will be
0. Guess a value for $f$, based on $\epsilon/D$ and a Reynolds number close to $\infty$.
1. Calculate $V$ from the energy equation with the current guess of $f$.
2. Calculate $Re_D = \rho V D/\mu$ from $V$.
3. Calculate a new value of $f$. Check if it is equal to the guess used in 1. If
yes, we stop. If not, then return to step 1 with this $f$ as our new guess.

Let's prepare ourselves for this iteration by defining some equations. The
roughness coefficient:
=#
eD = RoughnessCoefficient(ϵ/D)

# Here is $V$ as a function of $f$ (and the other quantities, which don't change).
Vfromf(f) = Velocity(sqrt(2*g*hp/(1+KL+f*L/D)))

# Here is $Re_D$ as a function of $V$ (and other quantities that don't change)
ReD(V) = ReynoldsNumber(ρ*V*D/μ)

#=
To do the iteration, we will define a simple function that takes in
an initial guess for $f$, iterates using the procedure above, and then returns
the correct value.

In the function, we create a quantity `f_old`, which simply hold our old guess
for $f$, and `f` will hold our next guess. We will compare these in each iteration,
checking whether they are nearlyequal. We'll say they're equal if they are closer
than $10^{-8}$. To ensure that the test fails at first, we initialize `f_old` to
infinity. This will force it to proceed into the iteration loop:
=#
function iterate_f(f0::FrictionFactor)
  f = f0
  f_old = Inf
  while abs(f - f_old) > 1e-8
      f_old = f
      f = FrictionFactor(ReD(Vfromf(f_old)),eD)
  end
  return f
end


#=
Now, to get $f$ from $Re_D$ and $\epsilon/D$, we use the function
`FrictionFactor(Re,eD)`, which solves the Colebrook equation for $f$.
Here, we use it to find our initial guess for $f$, treating $Re_D$ as very large.
=#
f0 = FrictionFactor(ReynoldsNumber(1e10),eD)

f = iterate_f(f0)

# We converged on the final $f$ value! The actual velocity is thus
V = Vfromf(f)
value(Vfromf(f),u"ft/s")

# And the final Reynolds number is
ReD(V)

# And finally, the volume flow rate after making the change to the setup is
Q = VolumeFlowRate(V*A)
value(Q,u"ft^3/s")

#=
So this shows that the flow rate is too small to meet the specifications.
We would need to decrease the length, $L$, or reduce the losses in $K_L$
to make it meet the specifications.
=#
