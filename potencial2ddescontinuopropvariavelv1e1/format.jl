function calcula_arco(x1,y1,x2,y2,xc,yc)
# Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
# horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
# and the horizontal direction


dx1 = x1 - xc; dy1 = y1 - yc
dx2 = x2 - xc; dy2 = y2 - yc

# Computation of tet1
      if dy1 == 0				# The point 1 and the center have the same y coordinate
        if x1 > xc
          tet1 = 0
        else  # (x1 < xc)
          tet1 = π
        end
      elseif dx1 == 0				# The point 1 and the center have the same x coordinate
        if y1 > yc
          tet1 = π/2
        else  # (y1 < yc)
          tet1 = -π/2
        end
      else  # (dx1~=0 e dy1~=0)
        tet1 = atan(dy1/dx1);
        if dx1<0 && tet1<0
          tet1 = π + tet1
        elseif dx1 < 0 && tet1>0
          tet1 = -π + tet1
        end
      end

# Computation of tet2
      if dy2 == 0				# The point 2 and the center have the same y coordinate
        if x2 > xc
          tet2 = 0
        else  # (x2 < xc)
          tet2 = π
        end
      elseif dx2 == 0				# The point 2 and the center have the same x coordinate
        if y2 > yc
          tet2 = π/2
        else  # (y2 < yc)
          tet2 = -π/2
        end
      else  # (dx2~=0 e dy2~=0)
        tet2 = atan(dy2/dx2);
        if dx2<0 && tet2<0
          tet2 = π + tet2
        elseif dx2 < 0 && tet2>0
          tet2 = -π + tet2
        end
      end
[tet1,tet2]
end


function calcula_centro(x1,y1,x2,y2,raio)
# Compute the center of an arc given two points and the radius

xm=(x1+x2)/2
ym=(y1+y2)/2
b=√((x2-x1)^2+(y2-y1)^2)
t1=(x2-x1)/b
t2=(y2-y1)/b
n1=t2
n2=-t1
h=√(abs(raio^2-(b/2)^2))
if(raio>0)
   if(n1==0)
      xc=xm
      yc=ym-n2/abs(n2)*h
   else
      xc=-n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
      yc=n2/n1*(xc-xm)+ym
   end
else
   if(n1==0)
      xc=xm
      yc=ym+n2/abs(n2)*h
   else
      xc=n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
      yc=n2/n1*(xc-xm)+ym
   end
end
[xc,yc]
end

function format_dad(PONTOS,SEGMENTOS,MALHA,tipo_elem,afasta,vizinhos=false)

  # Programa para formatação dos dados de entrada
  num_elementos = sum(MALHA[:,2])                             #soma todos elementos
  ELEM = zeros(UInt32,num_elementos,2+tipo_elem)
  NOS = zeros(0,2)

  cont_nos = UInt16(0);                      # Counter to the physical nodes
  cont_el = UInt16(0);	                     # Counter to the elements (the number of physical and geometric elements is the same).
  num_lin = length(SEGMENTOS[:,1]);	         # Número de segmentos no contorno
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
  while (t < num_lin)                                            # While para cada segmento
    while (p2 != p_ini)       #DÚVIDA:p2 nunca vai ser p_ini     # while para cada elemento dentro do while de cada segmento
      num_el_lin = MALHA[t,2];
      p1  = round(UInt64,SEGMENTOS[t,2])              # Coordinates of the initial and final PONTOS of each line
      p2  = round(UInt64,SEGMENTOS[t,3])
      x1l = PONTOS[p1,2]        # [x1l,y1l,x2l,y2l]
      y1l = PONTOS[p1,3]
      x2l = PONTOS[p2,2]
      y2l = PONTOS[p2,3]
      if (SEGMENTOS[t,4]==0)     # The segment is a straight line
        delta_x = x2l - x1l                                          # Increment in x and y direction
        delta_y = y2l - y1l
      else                       # The segment is an arc
        r = SEGMENTOS[t,4]
        xc,yc=calcula_centro(x1l,y1l,x2l,y2l,r)                      # Compute the center of the arc and its coordinates

        r1 = √((x1l-xc)^2+(y1l-yc)^2)                                # Distance between p1 and center (r1)
        r2 = √((x2l-xc)^2+(y2l-yc)^2)                                # and between p2 and center (r2)

        if abs(r1-r2) < 0.00001*max_dl
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
      qsi = range(0, stop = 1, length = tipo_elem)                    #gerando qsis (lista)
      for i = 1 : num_el_lin                                          #for p/ elementos no segmento
        if (SEGMENTOS[t,4] == 0) # The segment is a straight line
          x_i = x1l + delta_x/num_el_lin*(i-1)		# initial x coordinate of the element
          y_i = y1l + delta_y/num_el_lin*(i-1)		# initial y coordinate of the element
          x_f = x1l + delta_x/num_el_lin*(i)      # final x coordinate of the element
          y_f = y1l + delta_y/num_el_lin*(i)      # final y coordinate of the element
          xs= x_i.+(x_f-x_i)*qsi                                      #Nós do elemento (lista) em função de qsi  x(qsi)
          ys= y_i.+(y_f-y_i)*qsi                                      #Nós do elemento (lista) em função de qsi  y(qsi)
        else  # The segment is an arc
          # Compute the node coordinates
          xs = xc+r1*cos(tet1+(i-1 .+qsi)*sig*divtet)
          ys = yc+r1*sin(tet1+(i-1 .+qsi)*sig*divtet)

        end
        NOS = [NOS;[xs ys]]                                 #definindo matriz NÓS

        cont_el = cont_el + 1                               #começa em 0 
        nos = (cont_el-1)*tipo_elem+1:cont_el*tipo_elem     #tipo = 3(quadrat)
        ELEM[cont_el,:] = [cont_el nos' tipo_elem]          #definindo matriz ELEM [count el  nó1 nó2 nó3  tipo]
      end
      t=t+1
    end
  end
  if vizinhos                                           #DÚVIDA pra que serve vizinhos e pq recalcular NOS?
    NOSn = unique(NOS,dims=1)                           # elimina nós iguais
    ELEMn = zeros(size(ELEM,1),tipo_elem)
    for i = 1:size(NOSn,1)                              # for para cada nó distinto
      equalind = findall((NOSn[i,1].==NOS[:,1])+(NOSn[i,2].==NOS[:,2]).==2)   # true +true = 2
      for ind in equalind
        ELEMn[ELEM[:,2:(1+tipo_elem)].==ind].=i
      end
    end

    vizinho= zeros(UInt32,size(ELEM,1),2)          # matriz do elem ini e o elem fin nos quais estao sendo calculado o salto
 
    for el=1:size(ELEM,1)                          #for para cada elemento
     # @infiltrate                                #aqio calcula se uma matriz com os saltos
     vizinho[el,:] = [ELEM[ELEMn[el,1].==ELEMn[:,end],1] ELEM[ELEMn[el,end].==ELEMn[:,1],1]]; #busca o no adjacente ao no el e monta a matriz SALTOSI[el no_adjacente]
    end
  end
  
  for i=1:size(ELEM,1)
    qsis = range(-1+afasta,stop=1-afasta,length=tipo_elem)    # Parametrização de -1 a 1 com afasta
    N_geo = calc_fforma_gen.(qsis,tipo_elem,0)                #funções de forma generalizada
    xn=zeros(tipo_elem,2)
    for k=1:tipo_elem
      x = NOS[ELEM[i,2:end-1],:]           #coordenadas dos nos geometricos que compoem o elemento
      xn[k,:] = N_geo[k]*x                 #coordenadas geometricas do no em qsi = 0
    end
    NOS[ELEM[i,2:end-1],:]=xn
  end
  #_____
  if vizinhos
    return NOS,ELEM,vizinho
  else
    return  NOS,ELEM
  end
end


function monta_Teq(dad,x)
# Separa fluxo e temperatura

# ncdc = número de linhas da matriz CDC
# T = vetor que contêm as temperaturas nos nós
# q = vetor que contêm o fluxo nos nós

ncdc = length(dad.tipoCDC)
T=zeros(ncdc)
q=zeros(ncdc)
for i=1:ncdc # Laço sobre as condições de contorno
    if dad.tipoCDC[i] == 1 # Fluxo é conhecido
        T[i] = x[i]; # A temperatura é o valor calculado
        q[i] = dad.valorCDC[i]; # O fluxo é a condição de contorno
    else # A temperatura é conhecida
        T[i] = dad.valorCDC[i]; # A temperatura é a condição de contorno
        q[i] = x[i]; # O fluxo é o valor calculado
    end
end
return T,q
end
