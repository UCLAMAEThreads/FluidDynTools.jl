import ThermofluidQuantities: FrictionFactor

@nondimvar RoughnessCoefficient

function FrictionFactor(Re::ReynoldsNumber,eD::ET) where ET <: DimensionlessPhysicalQuantity
    xr = (5,Inf)
    xroot = find_zero(x -> -2.0*log10(eD/3.7 + 2.51/Re*x) - x,xr,order=16)
    FrictionFactor(1/xroot^2)
end
