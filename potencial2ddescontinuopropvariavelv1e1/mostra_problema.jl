function mostra_geometria(NOS_GEO,ELEM_GEO)
p=plot(leg=false);
for i=1:size(ELEM_GEO,1)
p=plot!(NOS_GEO[ELEM_GEO[i,2:end-1],1],NOS_GEO[ELEM_GEO[i,2:end-1],2],m=(:+,1),c=:black)
end
p
end


function mostra_geometria(dad)
p=plot(leg=false);
for i=1:size(dad.ELEM_GEO,1)
p=plot!(dad.NOS_GEO[dad.ELEM_GEO[i,2:end-1],1],dad.NOS_GEO[dad.ELEM_GEO[i,2:end-1],2],m=(:+,1),c=:black)
p=scatter!(dad.NOS[dad.ELEM[i,2:end-1],1],dad.NOS[dad.ELEM[i,2:end-1],2],marker = (:cross, 10,8.3, :red, stroke(:red)))
end
p
end


function mostra_resultado(dad,PONTOS_INT,T,Ti)
p=plot(leg=false);

pts=[dad.NOS;PONTOS_INT];
Ts=[T;Ti]
Tsr=(Ts.-minimum(Ts))/maximum(Ts)
tri=Triangle.basic_triangulation(pts,collect(1:size(pts,1)))
for i in tri
t=Shape([pts[i[1],1],pts[i[2],1],pts[i[3],1]],[pts[i[1],2],pts[i[2],2],pts[i[3],2]])
c=(Tsr[i[1]]+Tsr[i[2]]+Tsr[i[3]])/3
plot!(t,fill = (0, 0.5, RGB(c,c,c)))
end
p
end

function mostra_resultado_vtk(dad,PONTOS_INT,T,Ti,nome)

pts=[dad.NOS;PONTOS_INT];
Ts=[T;Ti]


tri=Triangle.basic_triangulation(pts,collect(1:size(pts,1)))
for i in tri
t=Shape([pts[i[1],1],pts[i[2],1],pts[i[3],1]],[pts[i[1],2],pts[i[2],2],pts[i[3],2]])
c=(Tsr[i[1]]+Tsr[i[2]]+Tsr[i[3]])/3
end

cells = MeshCell.(VTKCellTypes.VTK_TRIANGLE, tri)
vtkfile = vtk_grid(nome, points, cells)
vtkfile["my_point_data", VTKPointData()] = Ts
outfiles = vtk_save(vtkfile)

end
