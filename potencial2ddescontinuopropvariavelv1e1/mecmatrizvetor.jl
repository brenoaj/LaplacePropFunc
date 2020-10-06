function mecmatrizvetor(NOS,NOS_GEO,tipoCDC,valorCDC,normal,ELEM,kmat)

n=size(NOS,1); # Nï¿½mero de nï¿½s = Nï¿½mero de elementos

A=zeros(n,n); # Matriz A
b=zeros(1,n); # Vetor b
# NOS(2 NOS n)  = coordenadas dos extremos dos elementos
# (inï¿½cio e final dos elementos). NOS=[x1 do nï¿½ 1, x1 do nï¿½ 2, ...
                             #       x2 do nï¿½ 1; x2 do nï¿½ 2; ...
# NOS_GEO(2 NOS n) = coordenadas dos nï¿½s (pontos fontes)
# bc(2,n) = matriz de condiï¿½ï¿½es de contorno
# bc(1,i) = 0 => o deslocamento na direï¿½ï¿½o NOS do nï¿½ i ï¿½ conhecido
# bc(2,i) = 0 => o deslocamento na direï¿½ï¿½o NOS_GEO do nï¿½ i ï¿½ conhecido
# bc(1,i) = 1 => a forï¿½a de superfï¿½cie na direï¿½ï¿½o NOS do nï¿½ i ï¿½ conhecida
# bc(2,i) = 1 => a forï¿½a de superfï¿½cie na direï¿½ï¿½o NOS_GEO do nï¿½ i ï¿½ conhecida
ELEM=ELEM'
NOS=NOS'
NOS_GEO=NOS_GEO'
for j=1:n # Loop on elements (Column)
    al = √((NOS_GEO[1,ELEM[2,j]]-NOS_GEO[1,ELEM[1,j]])^2  +(NOS_GEO[2,ELEM[2,j]]-NOS_GEO[2,ELEM[1,j]])^2); # Comprimento do elemento
    for i=1:n # Loop on source points (Row)
        # Compute parameters used in the formulas for the two intergals
        x11 = NOS_GEO[1,ELEM[1,j]]-NOS[1,i]; # diferenï¿½a de coordenada NOS do
        # inï¿½cio do elemento atï¿½ o ponto fonte
        x21 = NOS_GEO[2,ELEM[1,j]]-NOS[2,i];# diferenï¿½a de coordenada NOS_GEO do
        # inï¿½cio do elemento atï¿½ o ponto fonte
        x12 = NOS_GEO[1,ELEM[2,j]]-NOS[1,i];# diferenï¿½a de coordenada NOS do
        # final do elemento atï¿½ o ponto fonte
        x22 = NOS_GEO[2,ELEM[2,j]]-NOS[2,i];# diferenï¿½a de coordenada NOS_GEO do
        # final do elemento atï¿½ o ponto fonte
        r1 = √(x11^2 + x21^2); # Distï¿½ncia do ponto fonte ao inï¿½cio
                       # do elemento
        r2 = √(x12^2 + x22^2); # Distï¿½ncia do ponto fonte ao final
                       # do elemento
        # normal = vetor normal n1 e n2
        d = x11*normal[1,j] + x21*normal[2,j]; # = delta x1 n1+delta y1 n2
                               # = dr1/dn = h
        t1 = -x11*normal[2,j] + x21*normal[1,j]; # =-delta x1 n2+delta y1 n1
                               # = dr1/dt
        t2 = -x12*normal[2,j] + x22*normal[1,j]; # =-delta x2 n2+delta y2 n1
                               # = dr1/dt
        ds = abs(d)
        dtheta = atan(ds*al,ds^2+t1*t2)
        # Elementos da matriz G
        g = (-dtheta*ds + al + t1*log(r1)-t2*log(r2))/(2π*kmat);

        if(d<0 )
            dtheta = -dtheta
        end
        # Elementos da matriz H
        h = dtheta/2π
        if(i==j)
            h = -0.5; # Diagonal da matriz H
        end
        if(tipoCDC[j]==0) # A temperatura ï¿½ conhecida
            A[i,j] = A[i,j]-g; # Os elementos de G vï¿½o para a matriz A
            b[i] = b[i] - h*valorCDC[j];# Os elementos de H vï¿½o para o vetor b
        else # O fluxo ï¿½ conhecido
            A[i,j] = A[i,j] + h; # Os elementos de H vï¿½o para a matriz A
            b[i] = b[i] + g*valorCDC[j];# Os elementos de G vï¿½o para o vetor b
        end
    end
end
return A,b
end

function domain_field(xfield,y,T,q,node,dnorm,kmat) #PONTOS_INT',NOS_GEO',T,q,ELEM',normal,k
# -------------------------------------------------------------------------------
#  Function used to evaluate the temperature/heat flux using the direct
#     evaluation of the integral representation
# -------------------------------------------------------------------------------
if (isempty(xfield))   # Pontos internos
  f=zeros(1)
  fx=zeros(1)
  fy=zeros(1)
else
    pi2 = π*2
    nfield=size(xfield,2)   # Número de pontos internos
    f=zeros(nfield)
    fx=zeros(nfield)
    fy=zeros(nfield)
    n=size(y,2)             # Número de pontos no contorno
    al= sqrt.((y[1,vec(node[2,:])]-y[1,vec(node[1,:])]).^2+ (y[2,vec(node[2,:])]-y[2,vec(node[1,:])]).^2)

    for j=1:n
        # Comprimento do elemento
        # Distâncias entre o ponto interno e os extremos dos elementos

        x11 = y[1,node[1,j]] .- xfield[1,:]
        x21 = y[2,node[1,j]] .- xfield[2,:]
        x12 = y[1,node[2,j]] .- xfield[1,:]
        x22 = y[2,node[2,j]] .- xfield[2,:]

        r1 =  sqrt.(x11.^2 + x21.^2); # Distância para o início do elemento
        r2 =  sqrt.(x12.^2 + x22.^2); # Distância para o final do elemento

        # Projeção do vetor distância no vetor normal ao elemento
        d  =  x11.*dnorm[1,j] + x21.*dnorm[2,j]; # Figura A.1, página 178
        t1 = -x11.*dnorm[2,j] + x21.*dnorm[1,j]; # Distância T1 da figura A.1
        t2 = -x12.*dnorm[2,j] + x22.*dnorm[1,j]; # Distância T2 da figura A.1
        ds = abs.(d)
        dtheta = atan.(ds.*al[j],ds.^2+t1.*t2)

        # Equação (A.5) do livro do Liu: elementos da matriz G
        g = @. -(-dtheta.*ds + al[j] + t1.*real(log.(r1))-t2.*real(log.(r2)))/(pi2*kmat)

        # Equação (A.7) com nx=1, ny=0, tx=0, ty=1.
        kkx = (dtheta.*dnorm[1,j] - real(log.(r2./r1)).*dnorm[2,j])/(pi2*kmat)
        # Equação(A.7) do livro do Liu com nx=0, ny=1, tx=-1, ty=0.
        kky = (dtheta.*dnorm[2,j] + real(log.(r2./r1)).*dnorm[1,j])/(pi2*kmat)
        hhx = @. -(-(t2./r2.^2-t1./r1.^2).*dnorm[1,j] - d.*(1.0/r2.^2-1.0/r1.^2).*dnorm[2,j])/pi2; # Equa��o (A.8) com nx=1
        # ny=0; tx=0; ty=1.
        hhy = @. -(-(t2./r2.^2-t1./r1.^2).*dnorm[2,j] + d.*(1.0/r2.^2-1.0/r1.^2).*dnorm[1,j])/pi2; # Equa��o (A.8) com nx=0
        #     if(d<=0)
        #         dtheta = -dtheta
        #     end
        dtheta = dtheta.*sign.(d)
        h = -dtheta/pi2; # Equação (A.6): Elementos da matriz
        f = f + g*q[j] - h*T[j]; # Integral (2.12) com o termo de
        # domínio igual a zero.
        fx = fx + kkx*q[j] - hhx*T[j]; # Integral (2.12) com o termo de
        # domínio igual a zero.
        fy = fy + kky*q[j] - hhy*T[j]; # Integral (2.12) com o termo de
        # domínio igual a zero.
    end
end
return f,fx*kmat,fy*kmat
end
