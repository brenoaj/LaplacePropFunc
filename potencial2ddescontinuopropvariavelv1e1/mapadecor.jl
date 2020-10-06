import GeometricalPredicates.getx
import GeometricalPredicates.gety
import GeometricalPredicates.getz
using  VoronoiDelaunay

type IndexedPoint2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _idx::Int64
    IndexedPoint2D(x, y, idx) = new(x, y, idx)
    IndexedPoint2D(x, y) = new(x, y, 0)
end

getx(p::IndexedPoint2D) = p._x
gety(p::IndexedPoint2D) = p._y
getidx(p::IndexedPoint2D) = p._idx

function mapadecor(shape,T,dT)
    npoints = size(shape, 1)
    tess = DelaunayTessellation2D{IndexedPoint2D}(npoints)
    min_value, max_value = (minimum(shape), maximum(shape))
    sc = (max_value - min_value) / (max_coord - min_coord)
    sc_shape = (shape .- min_value) / sc + min_coord
    points = IndexedPoint2D[IndexedPoint2D(sc_shape[i, 1], sc_shape[i, 2], i)
                            for i in 1:size(shape, 1)]
    push!(tess, points)
    t=collect(tess)[1:end]

    nt=length(t)
    index=zeros(nt,3)
  for i=1:nt
    index[i,1]=getidx(geta(t[i]))
    index[i,2]=getidx(getb(t[i]))
    index[i,3]=getidx(getc(t[i]))
  end

    #return index

    element_num= size(index,1)
    element_order=3
    cell_size = element_num * ( element_order + 1 )

    open("saida.vtk", "w") do f
      write(f, "# vtk DataFile Version 2.0\n")
      write(f, "saida.vtk\n")
      write(f, "ASCII\n")
      write(f, "\n")
      write(f, "DATASET UNSTRUCTURED_GRID\n")
      write(f, "POINTS $npoints double\n")
      writedlm(f,[shape zeros(npoints,1)])

      write(f, "\n" )
      write(f, "CELLS $element_num $cell_size \n")
      writedlm(f,[ ones(element_num,1)*element_order  index-1])

      write(f, "\n" )
      write(f, "CELL_TYPES $element_num \n")
      writedlm(f,ones(element_num,1)*5)

      write(f, "\n" )
      write(f, "POINT_DATA $npoints\n")
      write(f, "SCALARS Temperatura double\n" )
      write(f, "LOOKUP_TABLE default\n" )
      writedlm(f,T)
      write(f, "VECTORS dt double\n" )
      writedlm(f,[dT T*0])
    end
    return #
end
