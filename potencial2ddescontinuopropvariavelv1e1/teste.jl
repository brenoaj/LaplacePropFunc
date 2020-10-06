function format_dad(PONTOS,SEGMENTOS,MALHA,tipo_elem,afasta)

  # Programa para formata��o dos dados de entrada
  num_elementos=sum(MALHA[:,2])
  ELEM=zeros(UInt32,num_elementos,2+tipo_elem)
  NOS=zeros(0,2)

  cont_nos = UInt16(0);  # Counter to the physical nodes
  cont_el = UInt16(0);	# Counter to the elements (the number of physical and geometric elements is the same).
  num_lin = length(SEGMENTOS[:,1]);	# N�mero de linhas no contorno
  p_ini = round(UInt64,SEGMENTOS[1,2])

  #______________________________________________________________________
  # Definition of the biggest dimension of the problem
  max_dl = 0
  for lin = 1 : num_lin
    p1 = round(UInt64,SEGMENTOS[lin,2])
    p2 = round(UInt64,SEGMENTOS[lin,3])
    xp1 = PONTOS[p1,2]
    yp1 = PONTOS[p1,3]
    xp2 = PONTOS[p2,2]
    yp2 = PONTOS[p2,3]
    dl = √((xp1-xp2)^2+(yp1-yp2)^2)
    if dl > max_dl
      max_dl = dl
    end
  end
  #_____________________________________________________________________

  no_ini=1
  t=1
  p2=0
  no1_prox=0
  while(t<num_lin)  	# While over all lines
    while(p2!=p_ini)
      num_el_lin = MALHA[t,2];	# Number of the elements in the line t
      # Coordinates of the initial and final PONTOS of each line
      # [x1l,y1l,x2l,y2l)
      p1  = round(UInt64,SEGMENTOS[t,2])
      p2  = round(UInt64,SEGMENTOS[t,3])
      x1l = PONTOS[p1,2]
      y1l = PONTOS[p1,3]
      x2l = PONTOS[p2,2]
      y2l = PONTOS[p2,3]
      if(SEGMENTOS[t,4]==0) # The segment is a straight line
        # Increment in x and y direction
        delta_x = x2l - x1l
        delta_y = y2l - y1l
      else #The segment is an arc
        # Compute the center of the arc and its coordinates
        r = SEGMENTOS[t,4]
        xc,yc=calcula_centro(x1l,y1l,x2l,y2l,r)
        # Distance between p1 and c (r1) and between p2 and c (r2)
        r1 = √((x1l-xc)^2+(y1l-yc)^2)
        r2 = √((x2l-xc)^2+(y2l-yc)^2)
        if abs(r1-r2)<.00001*max_dl
          # Compute the angle between the lines from point c to p1 [tet1) and c to p2 (tet2]
          tet1,tet2 = calcula_arco(x1l,y1l,x2l,y2l,xc,yc)
          if tet2 < tet1
            tet2 = tet2 + 2*π
          end

          # Angle of the sector defined by the arc
          if SEGMENTOS[t,4] > 0
            tet = abs(tet2-tet1)
            sig = 1
          else
            tet = 2*π-abs(tet2-tet1)
            sig = -1
          end

          # Angle between two nodes of the line
          divtet = tet/(num_el_lin)
        else
          error("Error in the data input file: Wrong central point")
        end
      end
      # Generation of elements and nodes
      qsi=range(afasta,stop=1-afasta,length=tipo_elem)
      for i = 1 : num_el_lin

        if(SEGMENTOS[t,4]==0) # The segment is a straight line
          x_i = x1l + delta_x/num_el_lin*(i-1);			# initial x coordinate of the element
          y_i = y1l + delta_y/num_el_lin*(i-1);			# initial y coordinate of the element
          x_f = x1l + delta_x/num_el_lin*(i);	# midpoint x coordinate of the element
          y_f = y1l + delta_y/num_el_lin*(i);	# midpoint y coordinate of the element
          xs= x_i.+(x_f-x_i)*qsi
          ys= y_i.+(y_f-y_i)*qsi
        else  # The segment is an arc
          # Compute the node coordinates
          xs = xc+r1*cos(tet1+(i-1 .+qsi)*sig*divtet)
          ys = yc+r1*sin(tet1+(i-1 .+qsi)*sig*divtet)

        end
        NOS=[NOS;[xs ys]]

        cont_el = cont_el + 1
        nos=(cont_el-1)*tipo_elem+1:cont_el*tipo_elem
        ELEM[cont_el,:]=[cont_el nos' tipo_elem]
      end
      t=t+1
    end
  end
    return NOS,ELEM
  end
