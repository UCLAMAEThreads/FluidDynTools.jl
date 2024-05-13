### Boundary layer routines

function fsrhs!(du,u,p,t)
  _, m, _ = p

  du[1] = u[2]
  du[2] = u[3]
  du[3] = -0.5*(m+1)*u[1]*u[3] - m*(1-u[2]^2)
end

function fs_integrate(h0,p)
  Vw, m, ηmax = p

  f0 = -2Vw
  g0 = 0.0
  gL = 1.0

  u0 = [f0;g0;h0]
  ηspan = (0.0,ηmax)


  prob = ODEProblem(fsrhs!,u0,ηspan,p)
  sol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-15)

  resid = sol[2,end] - gL

  return resid, sol

end

struct FalknerSkan
  beta :: Float64
  m :: Float64
  eta :: Vector{Float64}
  psi :: Vector{Float64}
  u :: Vector{Float64}
  v :: Vector{Float64}
  omega :: Vector{Float64}
  d99 :: Float64
  dstar :: Float64
  theta :: Float64
  Cf :: Float64
end

"""
    falknerskan2(beta[,Vw=0][,etamax=10][,h0init=1.3]) -> FalknerSkan

Compute the Falkner-Skan boundary layer solution for parameter `beta`, for an external velocity \$U_e = A x^m\$,
where \$m = \\beta/(2-\\beta)\$

The parameter `beta` can be a value larger than -0.1988 (the separation case).
It returns a data structure of type `FalknerSkan` which contains the full solution
behavior, with the following fields

* `eta` -> the vertical coordinates \$\\eta = y/\\delta\$
* `u` -> the streamwise velocity profile (\$u/U_e\$)
* `psi` -> the streamfunction profile (\$\\psi/(U_e\\delta)\$)
* `v` -> the vertical velocity profile (\$v/U_e Re_x^{1/2}\$)
* `omega` -> the vorticity profile (\$\\omega_z \\delta/U_e\$)
* `d99` -> the proportionality factor on the 99 percent thickness (\$\\delta_{99}/\\delta\$)
* `dstar` -> the proportionality factor on the displacement thickness (\$\\delta^{*}/\\delta\$)
* `theta` -> the proportionality factor on the momentum thickness (\$\\theta/\\delta\$)
* `Cf` -> the proportionality factor on the skin friction coefficient (\$C_f Re_x^{1/2}\$)

where \$\\delta = \\sqrt{\\nu x/Ue} = x/Re_x^{1/2}\$.

For example,

```
julia> fs = falknerskan2(0);

julia> fs.Cf
0.6641146724303825
```
This indicates that the skin friction coefficient is approximately \$C_f = 0.664/Re_x^{1/2}\$.

If you want to specify a wall velocity, in dimensionless form \$V_w/U_e Re_x^{1/2}\$, then use the optional `Vw = ` keyword.

In some cases you may need to modify the initial guess for the \$f''(0)\$ for the shooting method. This is done
with the `h0init = ` keyword.

Similarly, in some cases you may wish to provide a larger maximum \$\\eta\$ coordinate, which you can do
with the `etamax = ` keyword.
"""
function falknerskan2(β;Vw = 0.0,etamax = 10.0,h0init=1.3)

  m = β/(2-β)
  p = [Vw,m,etamax]

  h0 = FluidDynTools.find_zero(x -> FluidDynTools.fs_integrate(x,p)[1],h0init)
  resid,sol = FluidDynTools.fs_integrate(h0,p)
  
  return FalknerSkan(β,m,sol.t,
                           sol[1,:],
                           sol[2,:],
                           _fs_v(sol,m),
                           -sol[3,:],
                           FluidDynTools.blthickness_99(sol),
                           FluidDynTools.blthickness_displacement(sol),
                           FluidDynTools.blthickness_momentum(sol),
                           FluidDynTools.blskinfriction(sol))
    
end

_fs_v(sol,m) = (1-m)/2*sol.t.*sol[2,:] - (1+m)/2*sol[1,:]

_fs_streamfunction(sol::OrdinaryDiffEq.ODESolution) = sol[1,:]
_fs_velocity(sol::OrdinaryDiffEq.ODESolution) = sol[2,:]
_fs_eta(sol::OrdinaryDiffEq.ODESolution) = sol.t

function blthickness_99(sol::OrdinaryDiffEq.ODESolution)

    u, η = _fs_velocity(sol), _fs_eta(sol)

    sign_diff = sign.(u .- 0.99)  # +1 where u > 0.99, -1 otherwise
    i0 = findfirst(diff(sign_diff) .== 2.0)
    return sol.t[i0]+ (0.99-u[i0])/(u[i0+1]-u[i0])*(η[i0+1]-η[i0])
end

function blthickness_displacement(sol::OrdinaryDiffEq.ODESolution)
    u, η = _fs_velocity(sol), _fs_eta(sol)

    uhalf = 0.5*(u[2:end].+u[1:end-1])
    return sum((1 .- uhalf).*diff(η))
end

function blthickness_momentum(sol::OrdinaryDiffEq.ODESolution)
    u, η = _fs_velocity(sol), _fs_eta(sol)

    uhalf = 0.5*(u[2:end].+u[1:end-1])
    return sum(uhalf.*(1 .- uhalf).*diff(η))
end

blskinfriction(sol::OrdinaryDiffEq.ODESolution) = 2*sol[3,1]
