using RecipesBase
using ColorTypes
using ViscousFlow
import PlotUtils: cgrad, palette, color_list
using LaTeXStrings

#=
@userplot Trajectories

@recipe function f(h::Trajectories)
  if length(h.args) > 2
      error("`trajectories` should be given one or two arguments.  Got: $(typeof(h.args))")
  end
  traj_array = h.args[1]

  xguide := L"x"
  yguide := L"y"
  aspect_ratio := 1
  size --> (700,400)

  if length(h.args) == 2
    sys = h.args[2]
    @series begin
      sys.base_cache.bl
    end
  end

  if isa(traj_array,ODESolution)
    @series begin
      linewidth := 2
      traj_array[1,:], traj_array[2,:]
    end
  else
    for traj in traj_array
      @series begin
        linewidth := 2
        traj[1,:], traj[2,:]
      end
    end
  end
end
=#

@userplot FieldTrajectory

@recipe function f(h::FieldTrajectory;fieldlabel="Field",deriv=0)
  if length(h.args) != 4
      error("`fieldtrajectory` should be given four arguments.  Got: $(typeof(h.args))")
  end
  traj, field, sys, pnum = h.args

  xtraj, ytraj = traj[pnum]

  layout := (2,1)
  size --> (600,600)
  xlims -> (-1,2)

  @series begin
    subplot := 1
    aspect_ratio := 1
    ylims --> (-0.5,0.5)
    title := "Particle trajectory"
    idxs := [pnum]
    traj
  end

  
  if length(sys) > 0
    body = surfaces(sys)[1]
    @series begin
      subplot := 1
      ylims --> (-0.5,0.5)
      xguide --> "x"
      yguide --> "y"
      aspect_ratio := 1
      surfaces(sys)
    end
    yb = -20
    yt = 20
    Xr = [minimum(body.x), maximum(body.x), maximum(body.x), minimum(body.x), minimum(body.x)]
    Yr = [yb,yb,yt,yt,yb]
    @series begin
      subplot := 2
      fillrange := 0
      fillalpha --> 0.2
      linealpha --> 0.2
      fillcolor := :lightgray
      linecolor := :lightgray
      label := "Body location"
      Xr,Yr
    end
  end

  
  if typeof(field) <: Tuple
    utraj,vtraj = field_along_trajectory(field,traj,pnum,deriv=deriv)
    minuv = min(minimum(utraj),minimum(vtraj)) - 0.5
    maxuv = max(maximum(utraj),maximum(vtraj)) + 0.5

    @series begin
      subplot := 2
      xguide := "x"
      label := "x $fieldlabel"
      xtraj, utraj
    end

    @series begin
      subplot := 2
      xguide := "x"
      label := "y $fieldlabel"
      title := "$fieldlabel components along trajectory"
      legend := true
      ylims := (minuv,maxuv)
      xtraj, vtraj
    end
  else
    straj = field_along_trajectory(field,traj,pnum,deriv=deriv)

    @series begin
      subplot := 2
      xguide := "x"
      label := "$fieldlabel"
      title := "$fieldlabel along trajectory"
      ylims := (minimum(straj)-0.5,maximum(straj)+0.5)
      legend := true
      xtraj, straj
    end
  end

end

@userplot Staticpressure

@recipe function f(h::Staticpressure;ambient=Pressure(1u"atm"),numarrows=20)
    if length(h.args) != 2
      error("`staticpressure` should be given two arguments.  Got: $(typeof(h.args))")
     end
    depths, fluids = h.args

    _depths = typeof(depths) <: AbstractArray ? value.(depths) : [value(depths)]
    _fluids = typeof(fluids) <: AbstractArray ? fluids : [fluids]
    p = Pressure[]
    d = Depth[]
    fs = Depth(0u"m")
    push!(p,ambient)
    push!(d,fs)
    prepend!(_depths,value(fs))
    dd = ustrip(_depths[end])/numarrows
    for i in eachindex(_fluids)
        di = Depth.(range(0,ustrip(_depths[i+1]-_depths[i]),step=dd)*unit(_depths[i+1]))
        pi = Pressure.(p[end] .+ SpecificWeight(_fluids[i])*di)
        append!(d,Depth.(d[end] .+ di))
        append!(p,pi)
    end

    flip := true
    xguide := "Pressure"
    yguide := "Depth"
    widen := false
    ymirror := true
    @series begin
        p, d
    end
    @series begin
        seriestype := Plots.quiver
        quiver := (-ustrip.(p),zeros(length(p)))
        p,d
    end
end
