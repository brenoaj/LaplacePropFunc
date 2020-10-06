function domain_T(PONTOS_INT,dad,T,q,npg=8,Pxy=PONTOS_INT) #NOS_GEO,ELEM_GEO,ELEM,kmat,afasta
    # Função para calcular temperatura interna
    n=size(PONTOS_INT,1)
    typ = eltype(Pxy)
    Ti = zeros(typ,n)
    nelem = size(dad.ELEM,1)
    qsi,w = gausslegendre(npg)    # Quadratura de gauss
    for i = 1:n
        pf = Pxy[i,:]    # coordenada (x,y) do ponto fonte (Real,Complex,Dual)
        pft = PONTOS_INT[i,:] # coordenada (x,y) do ponto fonte (Real)
        for j = 1:nelem
        	tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
            X = dad.NOS[dad.ELEM[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
            h = zeros(typ,1,tipo_elem)
            g = zeros(typ,1,tipo_elem)    # Array de 'tipo_elem' linhas

            eets=zeros(tipo_elem-1)
            bs=zeros(tipo_elem-1)
            qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem)
            for ii=1:tipo_elem-1
              rel=X[ii+1,:]-X[ii,:]
              eets[ii]=dot(rel,pft-X[ii,:])/norm(rel)^2 *(qsis[ii+1] - qsis[ii]) + qsis[ii]
              N=calc_fforma_gen(eets[ii],tipo_elem,dad.afasta)
              ps=N*X
              bs[ii]=norm(ps'- pft)
            end
            b,ind=findmin(bs)   # valor,índice
            eet=eets[ind]
            eta,Jt=sinhtrans(qsi,eet,b)     # Transformada
            # rel=X[end,:]-X[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
            # eet=2*dot(rel,pft-X[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
            # N_geo=calc_fforma_gen(eet,tipo_elem_geo,0)
            # ps=N_geo*X
            # b=norm(ps'-pft)
            # # b,ind=findmin(bs)
            # eta,Jt=sinhtrans(qsi,eet,b)
            for k = 1:npg
                N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
                dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
                pg = N*X    # Ponto de gauss interpolador
                r = pg'-pf      # Distancia entre

                dXdqsi = dN*X   # [dx/dξ dy/dξ]
                dgamadqsi = norm(dXdqsi)  # dΓ/dξ = J(ξ) Jacobiano
                sx = dXdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
                sy = dXdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
                nx = sy; # Componente x do vetor normal unitário
                ny =-sx; # Componente y do vetor normal unitário
                
                kponto=kf(pft)
                Qast,Tast = calsolfund(r,kponto,[nx,ny])
                h+= N*Qast*dgamadqsi*w[k]*Jt[k]*kponto
                g+= N*Tast*dgamadqsi*w[k]*Jt[k]*kponto
            end
            Ti[i] += (h*T[dad.ELEM[j,2:1+tipo_elem]] - g*q[dad.ELEM[j,2:1+tipo_elem]])[1]
        end
    end
    Ti
end

function dad1(PONTOS_INT) # Temperatura analitica
    typ = eltype(PONTOS_INT)
    n = size(PONTOS_INT,1)
    Ta = zeros(typ,n)
    for k = 1:n
        r = hypot(PONTOS_INT[k,1],PONTOS_INT[k,2])
        theta = acos(PONTOS_INT[k,1]/r)
        Ta[k] = sqrt(r)*cos(theta/2)
    end
    Ta
end

function dad1_d2T(PONTOS_INT)
    n = size(PONTOS_INT,1)
    dT = Array{Float64,2}(undef,size(PONTOS_INT,1),2)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    d2T= Array{Float64,2}(undef,size(PONTOS_INT,1),3)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    result = DiffResults.HessianResult(PONTOS_INT[1,:])  #ForwardDiff.DiffResults.HessianResult
    for i = 1:n
        ForwardDiff.hessian!(result,xy -> dad1(xy)[1],PONTOS_INT[i,:]')
        d2T[i,[1,3,2]] = reshape(DiffResults.hessian(result),(4,1))[[1,2,4]]
    end
    d2T
end

function d_AD(PONTOS_INT,dad,T,q,npg)
    Ti = Array{Float64,2}(undef,size(PONTOS_INT,1),1)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    dT = Array{Float64,2}(undef,size(PONTOS_INT,1),2)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    d2T= Array{Float64,2}(undef,size(PONTOS_INT,1),3)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    result = DiffResults.HessianResult(PONTOS_INT[1,:])  #ForwardDiff.DiffResults.HessianResult
    for i = 1:size(PONTOS_INT,1)
        ForwardDiff.hessian!(result,xy -> domain_T(PONTOS_INT[i,:]',dad,T,q,npg,xy)[1], [PONTOS_INT[i,1] PONTOS_INT[i,2]])
        Ti[i] = DiffResults.value(result)
        dT[i,:] = DiffResults.gradient(result)
        # d2T[i,1] = DiffResults.hessian(result)[1,1]
        # d2T[i,2] = DiffResults.hessian(result)[2,2]
        # d2T[i,3] = DiffResults.hessian(result)[1,2]
        d2T[i,[1,3,2]] = reshape(DiffResults.hessian(result),(4,1))[[1,2,4]]
    end

    return -dT*dad.k, d2T
end

function d_DF(PONTOS_INT,dad,T,q,npg,ϵ=1e-6) #(...,perturbação)
    # Gerando matrizes de saída
    dT = Array{Float64,2}(undef,size(PONTOS_INT))   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    d2T= Array{Float64,2}(undef,size(PONTOS_INT,1),3)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    d2T= Array{Float64,2}(undef,size(PONTOS_INT,1),3)   # d2T = zeros(size(PONTOS_INT),3)  #undef: memory garbage
    DL = maximum(dad.NOS_GEO, dims=1) - minimum(dad.NOS_GEO, dims=1) # maior segmento
    # Calculando com e sem perturbação
    Txy = domain_T(PONTOS_INT,dad,T,q,npg)      # Txy: Temperatura no ponto interno (x,y)
    Txmaish = domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .+ [ϵ*DL[1] 0]) #T(x+h,y)
    Txmenosh =domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .- [ϵ*DL[1] 0]) #T(x-h,y)
    Tymaish = domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .+ [0 ϵ*DL[2]]) #T(x,y+k)  h = k
    Tymenosh =domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .- [0 ϵ*DL[2]]) #T(x,y-k)
    Txymaish =domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .+ [ϵ*DL[1] ϵ*DL[2]]) #T(x,y-k)
    Txymenosh=domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .- [ϵ*DL[1] ϵ*DL[2]]) #T(x,y-k)
    # Primeira derivada
    dT[:,1] = (Txmaish - Txmenosh)/(2ϵ*DL[1])       #Implementação de diferenças finitas
    dT[:,2] = (Tymaish - Tymenosh)/(2ϵ*DL[2])
    # Segunda derivada
    d2T[:,1] = (Txmaish - 2Txy + Txmenosh)/(ϵ*DL[1])^2       #Implementação de diferenças finitas
    d2T[:,2] = (Tymaish - 2Txy + Tymenosh)/(ϵ*DL[2])^2  # df2 = (f(x+h) -2*f(x) + f(x-h))/h^2#
    d2T[:,3] = (Txymaish-Txmaish-Tymaish +2Txy -Txmenosh-Tymenosh+Txymenosh)/(2*DL[1]*DL[2]ϵ^2) # 2*e^2*prod(DL)

    return -dT*dad.k, d2T
end

function d_VC(PONTOS_INT,dad,T,q,npg,h = 1e-12)
    # Definindo variáveis para o calculo da derivada
    dT = Array{Float64,2}(undef,size(PONTOS_INT))
    d2T = Array{Float64,2}(undef,size(PONTOS_INT,1),3)
    DL = maximum(dad.NOS_GEO, dims=1) - minimum(dad.NOS_GEO, dims=1) # maior segmento
    #Perturbação
    Ti = domain_T(PONTOS_INT,dad,T,q,npg)
    Timx = domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .+[h*DL[1]*im 0])
    Timy = domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .+[0 h*DL[2]*im])
    Timxy = domain_T(PONTOS_INT,dad,T,q,npg,PONTOS_INT .+[h*DL[1]*im  h*DL[2]*im])
    # Primeira derivada
    dT[:,1] = imag(Timx)/(h*DL[1])
    dT[:,2] = imag(Timy)/(h*DL[2])
    # Segundo derivada
    d2T[:,1] = 2*(Ti - real(Timx))/(h*DL[1])^2  # vc2 = 2*(f(x) - real(f(x+im*h)) )/h^2
    d2T[:,2] = 2*(Ti - real(Timy))/(h*DL[2])^2
    d2T[:,3] = real(Timx + Timy - Ti -Timxy)/(DL[1]*DL[2]*h^2)
    return -dT*dad.k, d2T
end

function d_SF(PONTOS_INT,dad,T,q,npg=8)   # Calcula a primeira esegunda derivada no ponto interno
  n = size(PONTOS_INT,1)         # Quantidade de nos fisicos discretizando o contorno
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  nc = size(dad.NOS,1)
  Hx,Hy,Gx,Gy=zeros(n,nc),zeros(n,nc),zeros(n,nc),zeros(n,nc)
  Hxx,Hyy,Hxy,Gxx,Gyy,Gxy=zeros(n,nc),zeros(n,nc),zeros(n,nc),zeros(n,nc),zeros(n,nc),zeros(n,nc)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss

  for j=1:nelem  #Laço dos elementos
    tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
    x = dad.NOS[dad.ELEM[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos

    for i=1:n   #Laço dos pontos fontes
      pf= PONTOS_INT[i,:]  # coordenada
      hx,hy,gx,gy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)
      hxx,hyy,hxy,gxx,gyy,gxy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)
        # rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        # eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
        eets=zeros(tipo_elem-1)
        bs=zeros(tipo_elem-1)
        qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem)
        for ii=1:tipo_elem-1
          rel=x[ii+1,:]-x[ii,:]
          eets[ii]=dot(rel,pf-x[ii,:])/norm(rel)^2 *(qsis[ii+1] - qsis[ii]) + qsis[ii]
          N=calc_fforma_gen(eets[ii],tipo_elem,dad.afasta)
          ps=N*x
          bs[ii]=norm(ps'-pf)
        end
        b,ind=findmin(bs)
        eet=eets[ind]
        eta,Jt=sinhtrans(qsi,eet,b)

      # qsi,w = otimiza_npg(x,pf)
      # eta,Jt= telles(qsi,eet)
      for k = 1:npg
        N =calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
        pg =N*x    # Ponto de gauss interpolador
        r =pg'-pf      # Distancia entre
        dN =calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
        dxdqsi =dN*x   # dx/dξ & dy/dξ
        dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
        sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
        nx=sy; # Componente x do vetor normal unitário
        ny=-sx; # Componente y do vetor normal unitário
        dQdx,dTdx,dQdy,dTdy = calsolfund_dx_dy(r,dad.k,[nx,ny])
        d2Qdx,d2Tdx,d2Qdy,d2Tdy = calsolfund_d2x_d2y(r,dad.k,[nx,ny])
        dQdxdy,dTdxdy = calsolfund_dxdy(r,dad.k,[nx,ny])
        cont = N'*dgamadqsi*w[k]*Jt[k]
        hx +=cont*dQdx
        gx +=cont*dTdx
        hy +=cont*dQdy
        gy +=cont*dTdy
        hxx +=cont*d2Qdx
        gxx +=cont*d2Tdx
        hyy +=cont*d2Qdy
        gyy +=cont*d2Tdy
        hxy +=cont*dQdxdy
        gxy +=cont*dTdxdy
      end
      Hx[i,dad.ELEM[j,2:1+tipo_elem]] += hx
      Gx[i,dad.ELEM[j,2:1+tipo_elem]] += gx
      Hy[i,dad.ELEM[j,2:1+tipo_elem]] += hy
      Gy[i,dad.ELEM[j,2:1+tipo_elem]] += gy
      Hxx[i,dad.ELEM[j,2:1+tipo_elem]] += hxx
      Gxx[i,dad.ELEM[j,2:1+tipo_elem]] += gxx
      Hyy[i,dad.ELEM[j,2:1+tipo_elem]] += hyy
      Gyy[i,dad.ELEM[j,2:1+tipo_elem]] += gyy
      Hxy[i,dad.ELEM[j,2:1+tipo_elem]] += hxy
      Gxy[i,dad.ELEM[j,2:1+tipo_elem]] += gxy
    end
  end
  dTidx = Hx*T - Gx*q
  dTidy = Hy*T - Gy*q
  d2Tidx = Hxx*T - Gxx*q
  d2Tidy = Hyy*T - Gyy*q
  dTidxdy = Hxy*T - Gxy*q
  -dTidx*dad.k, -dTidy*dad.k,d2Tidx,d2Tidy,dTidxdy
end
function Qin_analitic(PONTOS_INT)        # função analitica
    npint = size(PONTOS_INT,1)
    Q_ANA = zeros(npint,2)
    for k = 1:npint
        r = hypot(PONTOS_INT[k,1],PONTOS_INT[k,2])
        theta = acos(PONTOS_INT[k,1]/r)
        Q_ANA[k,:] = -[cos(theta/2)/(2*sqrt(r))  sin(theta/2)/(2*sqrt(r))]
    end
    Q_ANA
end

# function d2_SF(PONTOS_INT,dad,T,q,npg=8)   # (NOS,NOS_GEO,ELEM,ELEM_GEO,kmat,afasta,npg=8)
#   n = size(PONTOS_INT,1)         # Quantidade de nos fisicos discretizando o contorno
#   nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
#   nc = size(dad.NOS,1)
#   Hx,Hy,Gx,Gy=zeros(n,nc),zeros(n,nc),zeros(n,nc),zeros(n,nc)
#   qsi,w = gausslegendre(npg)    # Quadratura de gauss
#
#   for j=1:nelem  #Laço dos elementos
#     tipo_elem_geo = dad.ELEM_GEO[j,end]   #Tipo do elemento geometrico
#     tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
#     x = dad.NOS_GEO[dad.ELEM_GEO[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
#
#     for i=1:n   #Laço dos pontos fontes
#       pf= PONTOS_INT[i,:]  # coordenada
#       hx,hy,gx,gy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)
#         # rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#         # eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
#         eets=zeros(tipo_elem-1)
#         bs=zeros(tipo_elem-1)
#         qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem)
#         for ii=1:tipo_elem-1
#           rel=x[ii+1,:]-x[ii,:]
#           eets[ii]=dot(rel,pf-x[ii,:])/norm(rel)^2 *(qsis[ii+1] - qsis[ii]) + qsis[ii]
#           N_geo=calc_fforma_gen(eets[ii],tipo_elem_geo,0)
#           ps=N_geo*x
#           bs[ii]=norm(ps'-pf)
#         end
#         b,ind=findmin(bs)
#         eet=eets[ind]
#         eta,Jt=sinhtrans(qsi,eet,b)
#
#       # qsi,w = otimiza_npg(x,pf)
#       # eta,Jt= telles(qsi,eet)
#       for k = 1:npg
#         N_geo =calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
#         N =calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
#         pg =N_geo*x    # Ponto de gauss interpolador
#         r =pg'-pf      # Distancia entre
#         dN_geo =calc_dfforma_gen(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
#         dxdqsi =dN_geo*x   # dx/dξ & dy/dξ
#         dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
#         sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
#         sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
#         nx=sy; # Componente x do vetor normal unitário
#         ny=-sx; # Componente y do vetor normal unitário
#         d2Qdx,d2Tdx,d2Qdy,d2Tdy = calsolfund_d2x_d2y(r,dad.k,[nx,ny])
#         hx +=N'*d2Qdx*dgamadqsi*w[k]*Jt[k]
#         gx +=N'*d2Tdx*dgamadqsi*w[k]*Jt[k]
#         hy +=N'*d2Qdy*dgamadqsi*w[k]*Jt[k]
#         gy +=N'*d2Tdy*dgamadqsi*w[k]*Jt[k]
#       end
#       Hx[i,dad.ELEM[j,2:1+tipo_elem]] += hx
#       Gx[i,dad.ELEM[j,2:1+tipo_elem]] += gx
#       Hy[i,dad.ELEM[j,2:1+tipo_elem]] += hy
#       Gy[i,dad.ELEM[j,2:1+tipo_elem]] += gy
#     end
#   end
#   d2Tidx = Hx*T - Gx*q
#   d2Tidy = Hy*T - Gy*q
#   d2Tidx, d2Tidy
# end

# function d_SF(PONTOS_INT,dad,T,q,npg=8)   # (NOS,NOS_GEO,ELEM,ELEM_GEO,kmat,afasta,npg=8)
#   n = size(PONTOS_INT,1)         # Quantidade de nos fisicos discretizando o contorno
#   nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
#   nc = size(dad.NOS,1)
#   Hx,Hy,Gx,Gy=zeros(n,nc),zeros(n,nc),zeros(n,nc),zeros(n,nc)
#   qsi,w = gausslegendre(npg)    # Quadratura de gauss
#
#   for j=1:nelem  #Laço dos elementos
#     tipo_elem_geo = dad.ELEM_GEO[j,end]   #Tipo do elemento geometrico
#     tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
#     x = dad.NOS_GEO[dad.ELEM_GEO[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
#
#     for i=1:n   #Laço dos pontos fontes
#       pf= PONTOS_INT[i,:]  # coordenada
#       hx,hy,gx,gy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)
#         # rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#         # eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
#         eets=zeros(tipo_elem-1)
#         bs=zeros(tipo_elem-1)
#         qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem)
#         for ii=1:tipo_elem-1
#           rel=x[ii+1,:]-x[ii,:]
#           eets[ii]=dot(rel,pf-x[ii,:])/norm(rel)^2 *(qsis[ii+1] - qsis[ii]) + qsis[ii]
#           N_geo=calc_fforma_gen(eets[ii],tipo_elem_geo,0)
#           ps=N_geo*x
#           bs[ii]=norm(ps'-pf)
#         end
#         b,ind=findmin(bs)
#         eet=eets[ind]
#         eta,Jt=sinhtrans(qsi,eet,b)
#
#       # qsi,w = otimiza_npg(x,pf)
#       # eta,Jt= telles(qsi,eet)
#       for k = 1:npg
#         N_geo =calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
#         N =calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
#         pg =N_geo*x    # Ponto de gauss interpolador
#         r =pg'-pf      # Distancia entre
#         dN_geo =calc_dfforma_gen(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
#         dxdqsi =dN_geo*x   # dx/dξ & dy/dξ
#         dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
#         sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
#         sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
#         nx=sy; # Componente x do vetor normal unitário
#         ny=-sx; # Componente y do vetor normal unitário
#         dQastdx,dTastdx,dQastdy,dTastdy = calsolfund_dx_dy(r,dad.k,[nx,ny])
#         hx +=N'*dQastdx*dgamadqsi*w[k]*Jt[k]
#         gx +=N'*dTastdx*dgamadqsi*w[k]*Jt[k]
#         hy +=N'*dQastdy*dgamadqsi*w[k]*Jt[k]
#         gy +=N'*dTastdy*dgamadqsi*w[k]*Jt[k]
#       end
#       Hx[i,dad.ELEM[j,2:1+tipo_elem]] += hx
#       Gx[i,dad.ELEM[j,2:1+tipo_elem]] += gx
#       Hy[i,dad.ELEM[j,2:1+tipo_elem]] += hy
#       Gy[i,dad.ELEM[j,2:1+tipo_elem]] += gy
#     end
#   end
#   dTidx = Hx*T - Gx*q
#   dTidy = Hy*T - Gy*q
#   dTidx, dTidy
# end


function Calc_HeG_interno(PONTOS_INT,dad,npg,tipo)#NOS_GEO,ELEM_GEO,ELEM,kmat,afasta
    # Função para calcular temperatura interna
    n=size(PONTOS_INT,1)
    nelem = size(dad.ELEM,1)
    qsi,w = gausslegendre(npg)    # Quadratura de gauss
    H=zeros(n,nelem*tipo)
    G=zeros(n,nelem*tipo)
    for i = 1:n
        pf = PONTOS_INT[i,:] # coordenada (x,y) do ponto fonte (Real)
        for j = 1:nelem
            tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
            X = dad.NOS[dad.ELEM[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
            h = zeros(1,tipo_elem)
            g = zeros(1,tipo_elem)    # Array de 'tipo_elem' linhas

            eets=zeros(tipo_elem-1)
            bs=zeros(tipo_elem-1)
            qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem)
            for ii=1:tipo_elem-1
              rel=X[ii+1,:]-X[ii,:]
              eets[ii]=dot(rel,pf-X[ii,:])/norm(rel)^2 *(qsis[ii+1] - qsis[ii]) + qsis[ii]
              N=calc_fforma_gen(eets[ii],tipo_elem,dad.afasta)
              ps=N*X
              bs[ii]=norm(ps'- pf)
            end
            b,ind=findmin(bs)   # valor,índice
            eet=eets[ind]
            eta,Jt=sinhtrans(qsi,eet,b)     # Transformada
            # rel=X[end,:]-X[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
            # eet=2*dot(rel,pft-X[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
            # N_geo=calc_fforma_gen(eet,tipo_elem_geo,0)
            # ps=N_geo*X
            # b=norm(ps'-pft)
            # # b,ind=findmin(bs)
            # eta,Jt=sinhtrans(qsi,eet,b)
            for k = 1:npg
                N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
                dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
                pg = N*X    # Ponto de gauss interpolador
                r = pg'-pf      # Distancia entre

                dXdqsi = dN*X   # [dx/dξ dy/dξ]
                dgamadqsi = norm(dXdqsi)  # dΓ/dξ = J(ξ) Jacobiano
                sx = dXdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
                sy = dXdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
                nx = sy; # Componente x do vetor normal unitário
                ny =-sx; # Componente y do vetor normal unitário

                kponto=kf(pf)                                                 #definir kno ponto fonte

                Qast,Tast = calsolfund(r,kponto,[nx,ny])
                # @infiltrate
                h+= N*Qast*dgamadqsi*w[k]*Jt[k]*kponto
                g+= N*Tast*dgamadqsi*w[k]*Jt[k]*kponto
            end
# @infiltrate
            H[i,dad.ELEM[j,2:1+tipo_elem]] += h'
            G[i,dad.ELEM[j,2:1+tipo_elem]] += g'
        end
    end
    H,G
end