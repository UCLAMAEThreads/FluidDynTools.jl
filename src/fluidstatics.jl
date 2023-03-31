## Fluid statics

export CircularPlate, SemicircularPlate, RectangularPlate, TriangularPlate

abstract type AbstractShape end
abstract type CircularPlate <: AbstractShape end
abstract type RectangularPlate <: AbstractShape end
abstract type SemicircularPlate <: AbstractShape end
abstract type TriangularPlate <: AbstractShape end


#ThermofluidQuantities.SpecificWeight(sg::SpecificGravity) = SpecificWeight(sg*SpecificWeight(Water))


##### Some geometric formulas #####
"""
    Area(::CircularPlate,radius)

Calculate the area of a circle, given the `radius`.
"""
ThermofluidQuantities.Area(::Type{CircularPlate},radius) = Area(π*radius^2)

"""
    Area(::SemicircularPlate,radius)

Calculate the area of a semicircle, given the `radius`.
"""
ThermofluidQuantities.Area(::Type{SemicircularPlate},radius) = Area(π*radius^2/2)

"""
    Area(::RectangularPlate,base,height)

Calculate the area of a rectangle, given the `base` and `height` lengths.
"""
ThermofluidQuantities.Area(::Type{RectangularPlate},base,height) = Area(base*height)

"""
    Area(::TriangularPlate,base,height,apex)

Calculate the area of a triangle, given the lengths `base`, `height`, and `apex`,
which describes the horizontal distance from the left vertex to the apex. The base is
assumed to lie along the ``x`` axis.
"""
ThermofluidQuantities.Area(::Type{TriangularPlate},base,height) = Area(base*height/2)


"""
    CentroidX(::CircularPlate,radius), CentroidY(::CircularPlate,radius)

Calculate the ``x`` and ``y`` coordinates of the centroid of a circle, given the `radius`.
The circle is assumed to be centered at the origin, so these simply return zero.
"""
CentroidX(::Type{CircularPlate},radius) = CentroidX(0)
CentroidY(::Type{CircularPlate},radius) = CentroidY(0)

"""
    CentroidX(::SemicircularPlate,radius), CentroidY(::SemicircularPlate,radius)

Calculate the ``x`` and ``y`` coordinates of the centroid of a semicircle, given the `radius`.
The semicircle is assumed to lie in the region above the ``x`` axis.
"""
CentroidX(::Type{SemicircularPlate},radius) = CentroidX(0)
CentroidY(::Type{SemicircularPlate},radius) = CentroidY(4*radius/3/π)

"""
    CentroidX(::RectangularPlate,base,height), CentroidY(::RectangularPlate,,base,height)

Calculate the ``x`` and ``y`` coordinates of the centroid of a rectangle, given the `base`
and `height`. The rectangle is assumed to be centered at the origin, so these simply return zero.
"""
CentroidX(::Type{RectangularPlate},base,height) = CentroidX(0)
CentroidY(::Type{RectangularPlate},base,height) = CentroidY(0)


"""
    CentroidX(::TriangularPlate,base,height,apex), CentroidY(::TriangularPlate,base,height,apex)

Calculate the area of a triangle, given the lengths `base`, `height`, and `apex`,
which describes the horizontal distance from the left vertex to the apex. The base is
assumed to lie along the ``x`` axis.
"""
CentroidX(::Type{TriangularPlate},base,height,apex) = CentroidX((base+apex)/3)
CentroidY(::Type{TriangularPlate},base,height,apex) = CentroidX(height/3)


"""
    SecondAreaMomentX(::CircularPlate,radius), SecondAreaMomentY(::CircularPlate,radius)


Calculate the second area moments of a circle about the centroid, given the `radius`.
"""
SecondAreaMomentX(::Type{CircularPlate},radius) = SecondAreaMomentX(π*radius^4/4)
SecondAreaMomentY(::Type{CircularPlate},radius) = SecondAreaMomentY(π*radius^4/4)

"""
    SecondAreaMomentX(::SemicircularPlate,radius), SecondAreaMomentY(::SemicircularPlate,radius)


Calculate the second area moments of a semicircle about the centroid, given the `radius`.
The semicircle is assumed to lie in the region above the ``x`` axis.
"""
SecondAreaMomentX(::Type{SemicircularPlate},radius) = SecondAreaMomentX((π/8-8/9/π)*radius^4)
SecondAreaMomentY(::Type{SemicircularPlate},radius) = SecondAreaMomentY(π*radius^4/8)

"""
    SecondAreaMomentX(::RectangularPlate,base,height), SecondAreaMomentY(::RectangularPlate,base,height)


Calculate the second area moments of a rectangle about the centroid, given the `base`
and `height` lengths. Note that ``x`` is parallel to the base and ``y`` to the height.
"""
SecondAreaMomentX(::Type{RectangularPlate},base,height) = SecondAreaMomentX(base*height^3/12)
SecondAreaMomentY(::Type{RectangularPlate},base,height) = SecondAreaMomentY(base^3*height/12)

"""
    SecondAreaMomentX(::TriangularPlate,base,height,apex), SecondAreaMomentY(::TriangularPlate,base,height,apex)

Calculate the second area moments of a triangle about the centroid, given the lengths `base`, `height`, and `apex`,
which describes the horizontal distance from the left vertex to the apex. The base is
assumed to lie along the ``x`` axis.
"""
SecondAreaMomentX(::Type{TriangularPlate},base,height,apex) = SecondAreaMomentX(base*height^3/36)
SecondAreaMomentY(::Type{TriangularPlate},base,height,apex) = SecondAreaMomentY((base^3*height-base^2*height*apex+base*height*apex^2)/36)

## Parallel axis theorem
"""
    SecondAreaMomentX(S,a...;x,y), SecondAreaMomentY(S,a...;x,y)

Calculate the second area moments about ``x`` and ``y``, using the
parallel axis theorem. If these arguments are omitted, it is assumed that
the moments are to be calculated about the centroid.
"""
SecondAreaMomentX(::Type{S},a...;x=CentroidX(S,a...),y=CentroidY(S,a...)) where S <: AbstractShape =
      SecondAreaMomentX(SecondAreaMomentX(S,a...) + Area(S,a...)*(CentroidY(S,a...)-y)^2)
SecondAreaMomentY(::Type{S},a...;x=CentroidX(S,a...),y=CentroidY(S,a...)) where S <: AbstractShape =
      SecondAreaMomentY(SecondAreaMomentY(S,a...) + Area(S,a...)*(CentroidX(S,a...)-x)^2)
