function cal_T(NOS,dad,T,q,npg=8,Pxy=dad.NOS,inos=0)
  n = size(NOS,1)         # Quantidade de nos fisicos discretizando o contorno
  typ = eltype(Pxy)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  H = zeros(typ,1,size(dad.NOS,1))
  G = zeros(typ,1,size(dad.NOS,1))

  for j = 1:nelem  #Laço dos elementos
     tipo_elem_geo = dad.ELEM_GEO[j,end]   #Tipo do elemento geometrico
     tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
     x = dad.NOS_GEO[dad.ELEM_GEO[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
     for i = 1:n   #Laço dos pontos fontes
       pf = Pxy[i,:]   # coordenada (x,y) do ponto fonte (Real,Complex,Dual)
       pft = NOS[i,:]  # coordenada
       h = zeros(typ,1,tipo_elem)
       g = zeros(typ,1,tipo_elem)    # Array de 'tipo_elem' linhas
       if n > 1 # não estiver executando um loop em AD ou VC
         inos = i
       end
       nosing = dad.ELEM[j,2:1+tipo_elem] .== inos   # Verifica se algum dos nos é ponto fonte
       if sum(nosing) == 1     # elemento contem o ponto fonte
         qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem) # Parametrização de -1 a 1
         eet=qsis[nosing]   #Determina a posição do pf no elemento  # do elemento com a distribuição dos nos
         eet=eet[1]     # transformando em float
         eta,Jt= telles(qsi,eet);
       else       # elemento Ñ contem o pf
         # eets=zeros(tipo_elem-1)
         # bs=zeros(tipo_elem-1)
         rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
         eet=2*dot(rel,pft-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
         N_geo=calc_fforma(eet,tipo_elem_geo,0)
         ps=N_geo*x
         b=norm(ps'-pft)
         # b,ind=findmin(bs)
         eta,Jt=sinhtrans(qsi,eet,b)
         # @show eta,Jt
       end
       # qsi,w = otimiza_npg(x,pf)
       # eta,Jt= telles(qsi,eet);
       for k = 1:npg
         # N_geo= calc_fforma(eta[k],tipo_elem_geo,0)  # funções de forma (geometrico)
         # N= calc_fforma(eta[k],tipo_elem,afasta)
         N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
         N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
         pg = N_geo*x    # Ponto de gauss interpolador
         r = pg'-pf      # Distancia entre

         # dN_geo = calc_dfforma(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
         dN_geo = calc_dfforma_gen(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
         dxdqsi = dN_geo*x   # dx/dξ & dy/dξ
         dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
         sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
         sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
         # @show sx
         Qast,Tast=calsolfund(r,dad.k,[sy,-sx])
         h+=N*Qast*dgamadqsi*w[k]*Jt[k]
         g+=N*Tast*dgamadqsi*w[k]*Jt[k]
       end
       H[1,dad.ELEM[j,2:1+tipo_elem]] = h'
       G[1,dad.ELEM[j,2:1+tipo_elem]] = g'
       # Tc[i] += 2*(h*T[ELEM[j,2:1+tipo_elem]] - g*q[ELEM[j,2:1+tipo_elem]])[1]
     end
   end
   H[inos]= 0.0
   H[inos]= 0.5-sum(H) # hipotese de corpo com temperatura constante
   # return H
   return 2*(H*T - G*q) # temperatura no contorno
end
# @time plot(CVc(dad,T,q,npg),leg = false)
# @time plot([ADc(dad,T,q,npg)],title="oi",leg = false)
# @time plot([ADc(dad,T,q,npg)],title="oi",leg = false)
# i = 2
# plot([SFc(dad,T,q,npg) #=CVc(dad,T,q,npg)[2:3:end,:]=# Qin_analitic(dad.NOS)])
# plot([CVc(dad,T,q,npg) Qin_analitic(dad.NOS)])

function ADc(NOS,dad,T,q,npg=8)
  n = size(NOS,1)
  Qc = zeros(n,2)
  for i = 1:n
    Qc[i,:] = ForwardDiff.gradient(xy -> cal_T(NOS[i,:]',dad,T,q,npg,xy,i)[1], [NOS[i,1] NOS[i,2]])
  end
  -Qc*dad.k
end

function CVc(NOS,dad,T,q,npg=8,h=1e-6)
  n = size(NOS,1)
  dT = Array{Float64,2}(undef,n,2)
  DL = maximum(dad.NOS_GEO,dims=1) - minimum(dad.NOS_GEO,dims=1)

  for i = 1:n
    Timx = cal_T(NOS[i,:]',dad,T,q,npg,NOS[i,:]' + [h*DL[1]*im 0],i)[1]
    Timy = cal_T(NOS[i,:]',dad,T,q,npg,NOS[i,:]' + [0 h*DL[2]*im],i)[1]
    dT[i,1] = imag(Timx)/(h*DL[1])
    dT[i,2] = imag(Timy)/(h*DL[2])
  end
  return dT*dad.k
end

function SFc(NOS,dad,T,q,npg=8)
  n = size(NOS,1)
  dT = Array{Float64,2}(undef,n,2)
  for i = 1:n
    dT[i,:] = cal_TSF(NOS[i,:]',dad,T,q,npg,i)
  end
  -dT*dad.k
end


function cal_TSF(NOS,dad,T,q,npg=8,inos=0)
  n = size(NOS,1)         # Quantidade de nos fisicos discretizando o contorno
  typ = eltype(NOS)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  Hx,Hy = zeros(typ,1,size(dad.NOS,1)),zeros(typ,1,size(dad.NOS,1))
  Gx,Gy = zeros(typ,1,size(dad.NOS,1)),zeros(typ,1,size(dad.NOS,1))

  for j = 1:nelem  #Laço dos elementos
     tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
     x = dad.NOS_GEO[dad.ELEM[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
     for i = 1:n   #Laço dos pontos fontes
       pf = NOS[i,:]  # coordenada
       hx,hy = zeros(typ,1,tipo_elem),zeros(typ,1,tipo_elem)
       gx,gy = zeros(typ,1,tipo_elem),zeros(typ,1,tipo_elem)    # Array de 'tipo_elem' linhas
       if n > 1 # não estiver executando um loop em AD ou VC
         inos = i
       end
       nosing = dad.ELEM[j,2:1+tipo_elem] .== inos   # Verifica se algum dos nos é ponto fonte
       if sum(nosing) == 1     # elemento contem o ponto fonte
         qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem) # Parametrização de -1 a 1
         eet=qsis[nosing]   #Determina a posição do pf no elemento  # do elemento com a distribuição dos nos
         eet=eet[1]     # transformando em float
         eta,Jt= telles(qsi,eet);
       else       # elemento Ñ contem o pf
         # eets=zeros(tipo_elem-1)
         # bs=zeros(tipo_elem-1)
         rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
         eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
         N=calc_fforma(eet,tipo_elem,0)
         ps=N*x
         b=norm(ps'-pf)
         # b,ind=findmin(bs)
         eta,Jt=sinhtrans(qsi,eet,b)
         # @show eets
         # @show eta,Jt
       end
       # qsi,w = otimiza_npg(x,pf)
       # eta,Jt= telles(qsi,eet);
       for k = 1:npg
         # N_geo= calc_fforma(eta[k],tipo_elem_geo,0)  # funções de forma (geometrico)
         # N= calc_fforma(eta[k],tipo_elem,afasta)
        #  N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
         N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
         pg = N*x    # Ponto de gauss interpolador
         r = pg'-pf      # Distancia entre

         # dN_geo = calc_dfforma(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
         dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
         dxdqsi = dN*x   # dx/dξ & dy/dξ
         dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
         sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
         sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
         # @show sx
         dQastdx,dTastdx,dQastdy,dTastdy =calsolfund_dx_dy(r,dad.k,[sy,-sx])
         hx +=N*dQastdx*dgamadqsi*w[k]*Jt[k]
         gx +=N*dTastdx*dgamadqsi*w[k]*Jt[k]
         hy +=N*dQastdy*dgamadqsi*w[k]*Jt[k]
         gy +=N*dTastdy*dgamadqsi*w[k]*Jt[k]
       end
       Hx[1,dad.ELEM[j,2:1+tipo_elem]] = hx'
       Gx[1,dad.ELEM[j,2:1+tipo_elem]] = gx'
       Hy[1,dad.ELEM[j,2:1+tipo_elem]] = hy'
       Gy[1,dad.ELEM[j,2:1+tipo_elem]] = gy'
       # Tc[i] += 2*(h*T[ELEM[j,2:1+tipo_elem]] - g*q[ELEM[j,2:1+tipo_elem]])[1]
     end
   end
   Hx[inos],Hy[inos]= 0.0,0.0
   Hx[inos],Hy[inos]= -sum(Hx),-sum(Hy) # hipotese de corpo com temperatura constante

   return [2*(Hx*T - Gx*q)[1], 2*(Hy*T - Gy*q)[1]]
end


function cal_dTSF(NOS,dad,T,q,npg=32)
  n = size(NOS,1)         # Quantidade de nos fisicos discretizando o contorno
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  Hx,Hy = zeros(size(dad.NOS,1)),zeros(size(dad.NOS,1))
  Gx,Gy = zeros(size(dad.NOS,1)),zeros(size(dad.NOS,1))
  dt=zeros(size(NOS,1),2)
  for i = 1:n   #Laço dos pontos fontes
      pf = NOS[i,:]  # coordenada
       
      for j = 1:nelem  #Laço dos elementos

          tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
          hx,hy = zeros(tipo_elem),zeros(tipo_elem)
          gx,gy = zeros(tipo_elem),zeros(tipo_elem)    # Array de 'tipo_elem' linhas
    
          x = dad.NOS[dad.ELEM[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
          nosing = dad.ELEM[j,2:1+tipo_elem] .== i   # Verifica se algum dos nos é ponto fonte
          if sum(nosing) == 1     # elemento contem o ponto fonte
              qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem) # Parametrização de -1 a 1
              eet=qsis[nosing]   #Determina a posição do pf no elemento  # do elemento com a distribuição dos nos
              eet=eet[1]     # transformando em float
              #   eta,Jt= telles(qsi,eet);
              eta = cos.(pi*(0:2*npg+1)/(2*npg+1));
              hxv= zeros(2*npg+2,tipo_elem)  
              gxv= zeros(2*npg+2,tipo_elem)  
              hyv= zeros(2*npg+2,tipo_elem)  
              gyv= zeros(2*npg+2,tipo_elem)
              
              hxv1= zeros(2*npg+2,tipo_elem)
              hyv1= zeros(2*npg+2,tipo_elem)
              gxv1= zeros(2*npg+2,tipo_elem)
              gyv1= zeros(2*npg+2,tipo_elem)
              for k = 1:2*npg+2
                  # N_geo= calc_fforma(eta[k],tipo_elem_geo,0)  # funções de forma (geometrico)
                  # N= calc_fforma(eta[k],tipo_elem,afasta)
                  # N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
                  N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
                  pg = N*x    # Ponto de gauss interpolador
                  r = pg'-pf      # Distancia entre
                  
                  # dN_geo = calc_dfforma(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
                  dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
                  dxdqsi = dN*x   # dx/dξ & dy/dξ
                  dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
                  sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
                  sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
                  # @show sx
                  dQastdx,dTastdx,dQastdy,dTastdy =calsolfund_dx_dy(r,dad.k,[sy,-sx])
                  # @show N*dQastdx*dgamadqsi,hxv[k,:]
                                                                              
                  hxv[k,:] =N*dQastdx*dgamadqsi*(eta[k]-eet)^2
                  hyv[k,:] =N*dQastdy*dgamadqsi*(eta[k]-eet)^2
                  gxv[k,:] =N*dTastdx*dgamadqsi*(eta[k]-eet)^1
                  gyv[k,:] =N*dTastdy*dgamadqsi*(eta[k]-eet)^1
              end
              

              hx =clenshaw_curtis2(hxv,eet)
              hy =clenshaw_curtis2(hyv,eet)
              gx =clenshaw_curtis1(gxv,eet)
              gy =clenshaw_curtis1(gyv,eet)
          else       # elemento Ñ contem o pf
              # eets=zeros(tipo_elem-1)
              # bs=zeros(tipo_elem-1)
              rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
              eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
              N=calc_fforma(eet,tipo_elem,dad.afasta)
              ps=N*x
             b=norm(ps'-pf)
          #     # b,ind=findmin(bs)
              # eta1,Jt1=sinhtrans(qsi,eet,b)
              # eta,Jt=sinhtrans(eta1,eet,b)
              # Jt=Jt.*Jt1
               eta,Jt=telles(qsi,eet)

              # @show eets
              # @show eta,Jt
              for k = 1:npg
                  # N_geo= calc_fforma(eta[k],tipo_elem_geo,0)  # funções de forma (geometrico)
                  # N= calc_fforma(eta[k],tipo_elem,afasta)
                  # N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
                  N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
                  pg = N*x    # Ponto de gauss interpolador
                  r = pg'-pf      # Distancia entre
                  
                  # dN_geo = calc_dfforma(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
                  dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
                  dxdqsi = dN*x   # dx/dξ & dy/dξ
                  dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
                  sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
                  sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
                  # @show sx
                  dQastdx,dTastdx,dQastdy,dTastdy =calsolfund_dx_dy(r,dad.k,[sy,-sx])
                  hx +=N'*dQastdx*dgamadqsi*w[k]*Jt[k]
                  gx +=N'*dTastdx*dgamadqsi*w[k]*Jt[k]
                  hy +=N'*dQastdy*dgamadqsi*w[k]*Jt[k]
                  gy +=N'*dTastdy*dgamadqsi*w[k]*Jt[k]
              end
          end
          
          
          Hx[dad.ELEM[j,2:1+tipo_elem]] = hx'
          Gx[dad.ELEM[j,2:1+tipo_elem]] = gx'
          Hy[dad.ELEM[j,2:1+tipo_elem]] = hy'
          Gy[dad.ELEM[j,2:1+tipo_elem]] = gy'
          # Tc[i] += 2*(h*T[ELEM[j,2:1+tipo_elem]] - g*q[ELEM[j,2:1+tipo_elem]])[1]
      end
      # @show dot(Hx,T),dot(Gx,q),dot(Hy,T),dot(Gy,q)
      # @show Hx
      # Hx[i],Hy[i]= 0.0,0.0
      # Hx[i],Hy[i]= -sum(Hx),-sum(Hy) # hipotese de corpo com temperatura constante
      dt[i,:]= [2*(dot(Hx,T) - dot(Gx,q)), 2*(dot(Hy,T) - dot(Gy,q))]
  end
  #  
  #  
  
  return dt
end
function cal_d2TSF(NOS,dad,T,q,npg=32)
  n = size(NOS,1)         # Quantidade de nos fisicos discretizando o contorno
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  Hxx,Hyy,Hxy = zeros(size(dad.NOS,1)),zeros(size(dad.NOS,1)),zeros(size(dad.NOS,1))
  Gxx,Gyy,Gxy = zeros(size(dad.NOS,1)),zeros(size(dad.NOS,1)),zeros(size(dad.NOS,1))
  dt=zeros(size(NOS,1),3)
  for i = 1:n   #Laço dos pontos fontes
      pf = NOS[i,:]  # coordenada
       
      for j = 1:nelem  #Laço dos elementos

          # tipo_elem_geo = dad.ELEM_GEO[j,end]   #Tipo do elemento geometrico
          tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
          hxx,hyy,hxy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)
          gxx,gyy,gxy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)    # Array de 'tipo_elem' linhas
    
          x = dad.NOS[dad.ELEM[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
          nosing = dad.ELEM[j,2:1+tipo_elem] .== i   # Verifica se algum dos nos é ponto fonte
          if sum(nosing) == 1     # elemento contem o ponto fonte
              qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem) # Parametrização de -1 a 1
              eet=qsis[nosing]   #Determina a posição do pf no elemento  # do elemento com a distribuição dos nos
              eet=eet[1]     # transformando em float
              #   eta,Jt= telles(qsi,eet);
              eta = cos.(pi*(0:2*npg+1)/(2*npg+1));
              hxxv,hyyv,hxyv = zeros(2*npg+2,tipo_elem),zeros(2*npg+2,tipo_elem),zeros(2*npg+2,tipo_elem) 
              gxxv,gyyv,gxyv = zeros(2*npg+2,tipo_elem),zeros(2*npg+2,tipo_elem),zeros(2*npg+2,tipo_elem) 

              for k = 1:2*npg+2
                  # N_geo= calc_fforma(eta[k],tipo_elem_geo,0)  # funções de forma (geometrico)
                  # N= calc_fforma(eta[k],tipo_elem,afasta)
                  # N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
                  N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
                  pg = N*x    # Ponto de gauss interpolador
                  r = pg'-pf      # Distancia entre
                  
                  # dN_geo = calc_dfforma(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
                  dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta)# calcula dN\dξ N1,N2 e N3
                  dxdqsi = dN*x   # dx/dξ & dy/dξ
                  dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
                  sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
                  sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
                  # @show sx
                  d2Qdx,d2Tdx,d2Qdy,d2Tdy = calsolfund_d2x_d2y(r,dad.k,[sy,-sx])
                  dQdxdy,dTdxdy = calsolfund_dxdy(r,dad.k,[sy,-sx])

                  # @show N*dQastdx*dgamadqsi,hxv[k,:]
                  hxxv[k,:] =N*d2Qdx*dgamadqsi*(eta[k]-eet)^3
                  gxxv[k,:] =N*d2Tdx*dgamadqsi*(eta[k]-eet)^2
                  hyyv[k,:] =N*d2Qdy*dgamadqsi*(eta[k]-eet)^3
                  gyyv[k,:] =N*d2Tdy*dgamadqsi*(eta[k]-eet)^2
                  hxyv[k,:] =N*dQdxdy*dgamadqsi*(eta[k]-eet)^3
                  gxyv[k,:] =N*dTdxdy*dgamadqsi*(eta[k]-eet)^2
              end
              hxx =clenshaw_curtis3(hxxv,eet)
              gxx =clenshaw_curtis2(gxxv,eet)
              hyy =clenshaw_curtis3(hyyv,eet)
              gyy =clenshaw_curtis2(gyyv,eet)
              hxy =clenshaw_curtis3(hxyv,eet)
              gxy =clenshaw_curtis2(gxyv,eet)
          else       # elemento Ñ contem o pf
              # eets=zeros(tipo_elem-1)
              # bs=zeros(tipo_elem-1)
              rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
              eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
              N=calc_fforma(eet,tipo_elem,dad.afasta)
              ps=N*x
             b=norm(ps'-pf)
          #     # b,ind=findmin(bs)
              eta1,Jt1=sinhtrans(qsi,eet,b)
              eta,Jt=sinhtrans(eta1,eet,b)
              Jt=Jt.*Jt1
              #  eta,Jt=telles(qsi,eet)

              # @show eets
              # @show eta,Jt
              for k = 1:npg
                  # N_geo= calc_fforma(eta[k],tipo_elem_geo,0)  # funções de forma (geometrico)
                  # N= calc_fforma(eta[k],tipo_elem,afasta)
                  # N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0) #funções de forma generalizada
                  N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)
                  pg = N*x    # Ponto de gauss interpolador
                  r = pg'-pf      # Distancia entre
                  
                  # dN_geo = calc_dfforma(eta[k],tipo_elem_geo,0) # calcula dN\dξ N1,N2 e N3
                  dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta) # calcula dN\dξ N1,N2 e N3
                  dxdqsi = dN*x   # dx/dξ & dy/dξ
                  dgamadqsi =norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
                  sx=dxdqsi[1]/dgamadqsi; # Componente x do vetor tangente dx/dΓ
                  sy=dxdqsi[2]/dgamadqsi; # Componente y do vetor tangente dy/dΓ
                  # @show sx
                  d2Qdx,d2Tdx,d2Qdy,d2Tdy = calsolfund_d2x_d2y(r,dad.k,[sy,-sx])
                  dQdxdy,dTdxdy = calsolfund_dxdy(r,dad.k,[sy,-sx])                    

                  hxx +=N'*d2Qdx*dgamadqsi*w[k]*Jt[k]
                  gxx +=N'*d2Tdx*dgamadqsi*w[k]*Jt[k]
                  hyy +=N'*d2Qdy*dgamadqsi*w[k]*Jt[k]
                  gyy +=N'*d2Tdy*dgamadqsi*w[k]*Jt[k]
                  hxy +=N'*dQdxdy*dgamadqsi*w[k]*Jt[k]
                  gxy +=N'*dTdxdy*dgamadqsi*w[k]*Jt[k]
              end
          end
          
          
          Hxx[dad.ELEM[j,2:1+tipo_elem]] += hxx
          Gxx[dad.ELEM[j,2:1+tipo_elem]] += gxx
          Hyy[dad.ELEM[j,2:1+tipo_elem]] += hyy
          Gyy[dad.ELEM[j,2:1+tipo_elem]] += gyy
          Hxy[dad.ELEM[j,2:1+tipo_elem]] += hxy
          Gxy[dad.ELEM[j,2:1+tipo_elem]] += gxy
          # Tc[i] += 2*(h*T[ELEM[j,2:1+tipo_elem]] - g*q[ELEM[j,2:1+tipo_elem]])[1]
      end
      #  @show dot(Hxy,T),dot(Gxy,q),dot(Hxx,T),dot(Gxx,q),dot(Hyy,T),dot(Gyy,q)
      # @show Hxy
      # @show i,sum(Hxx),sum(Hyy),sum(Hxy)
      # @show Gxy
      dt[i,:]= [2*(dot(Hxx,T) - dot(Gxx,q)), 2*(dot(Hyy,T) - dot(Gyy,q)), 2*(dot(Hxy,T) - dot(Gxy,q))]
  end
  #  Hx[inos],Hy[inos]= 0.0,0.0
  #  Hx[inos],Hy[inos]= -sum(Hx),-sum(Hy) # hipotese de corpo com temperatura constante
  
  return dt
end


function clenshaw_curtis(f)         
  n=size(f,1)-1
  fx = f/(2*n); 
  g = real(fft(fx[[1:n+1;n:-1:2]]));
  a = [g[1]; g[2:n]+g[2*n:-1:n+2]; g[n+1]];
  w = 0*a; w[1:2:end] = 2 ./(1 .-(0:2:n).^2);
  I = w'*a;
end

function clenshaw_curtis1(f,e)         
  n=size(f,1)-1
  I=zeros(size(f,2))
  for i=1:size(f,2)
  fx = f[:,i]/(2*n);     g = real(fft(fx[[1:n+1;n:-1:2]]));
  a = [g[1]; g[2:n]+g[2*n:-1:n+2]; g[n+1]];
  d= 0*a
  for k=0:n-1
      if k==0
          d[n-k]=a[n+1-k]+2*e*d[n+1-k]
      else
          d[n-k]=2*a[n+1-k]+2*e*d[n+1-k]-d[n+2-k]   
      end   
  end
  w = 0*a; w[1:2:end] = 2 ./(1 .-(0:2:n).^2);w[1]=1;
  fe=a'cos.((0:n)*acos(e))
  I[i] = w'*d+fe*log((1-e)/(1+e));
end
I
end

function clenshaw_curtis2(f,e)         
  n=size(f,1)-1
  I=zeros(size(f,2))
  for i=1:size(f,2)
  fx = f[:,i]/(2*n); 
  g = real(fft(fx[[1:n+1;n:-1:2]]));
  a = [g[1]; g[2:n]+g[2*n:-1:n+2]; g[n+1]];
  d= 0*a
  dl= 0*a
  for k=0:n-1
      if k==0
          d[n-k]=a[n+1-k]+2*e*d[n+1-k]
          dl[n-k]=d[n+1-k]+2*e*dl[n+1-k]
      else
          d[n-k]=2*a[n+1-k]+2*e*d[n+1-k]-d[n+2-k]   
          dl[n-k]=2*d[n+1-k]+2*e*dl[n+1-k]-dl[n+2-k]
      end   
  end
  w = 0*a; w[1:2:end] = 2 ./(1 .-(0:2:n).^2);w[1]=1;
  fle=a'*((0:n).* sin.((0:n)*acos(e)))/sqrt(1 - e^2)
  fe=a'cos.((0:n)*acos(e))
  I[i] = w'*dl+fle*log((1-e)/(1+e))-2*fe/(1-e^2);
end
I
end

function clenshaw_curtis3(f,e)         
  n=size(f,1)-1
  I=zeros(size(f,2))
  for i=1:size(f,2)
  fx = f[:,i]/(2*n); 
  g = real(fft(fx[[1:n+1;n:-1:2]]));
  a = [g[1]; g[2:n]+g[2*n:-1:n+2]; g[n+1]];
  d= 0*a
  dl= 0*a
  dl2= 0*a
  for k=0:n-1
      if k==0
          d[n-k]=a[n+1-k]+2*e*d[n+1-k]
          dl[n-k]=d[n+1-k]+2*e*dl[n+1-k]
          dl2[n-k]=dl[n+1-k]+2*e*dl2[n+1-k]
      else
          d[n-k]=2*a[n+1-k]+2*e*d[n+1-k]-d[n+2-k]   
          dl[n-k]=2*d[n+1-k]+2*e*dl[n+1-k]-dl[n+2-k]
          dl2[n-k]=2*dl[n+1-k]+2*e*dl2[n+1-k]-dl2[n+2-k]
      end   
  end
  w = 0*a; w[1:2:end] = 2 ./(1 .-(0:2:n).^2);w[1]=1;
  # @show n,e
  f2le=a'*((collect(0:n).*e.*sin.((0:n).*acos(e)))/(1-e^2)^(3/2)-((0:n).^2 .*cos.((0:n).*acos(e)))/(1-e^2))
  fle=a'*((0:n).* sin.((0:n)*acos(e)))/sqrt(1 - e^2)
  fe=a'cos.((0:n)*acos(e))
  I[i] = (w'*dl+f2le*log((1-e)/(1+e))-4*fle/(1-e^2)-4*e*fe/(1-e^2)^2)/2+f2le/8;
end
I
end


"Funcao para suavizar"
function suavizar(dad,vs,T)
  ndelinhas_ELEM = size(dad.ELEM,1) #número de elementos geometricos na malha
  Ts = zeros(Float64,ndelinhas_ELEM,2) #cria a matriz para armazenar a temperatura nos extremos do elemento. Ts = [T(qsi = -1), T(qsi = 1)].
Tsuave=0*T
  #Calculo dos saltos nos elementos
  for el = 1 : ndelinhas_ELEM
    tipo = dad.ELEM[el,end] #define o tipo do elemento
    Ne = calc_fforma_gen(-1.,tipo,dad.afasta); #funcao de forma para o no adjacente a esquerda
    Nd = calc_fforma_gen(1.,tipo,dad.afasta); #funcao de forma para o no adjacente a direita
    Tn = T[dad.ELEM[el,2:end-1]]; #busca o valor nodal da temperatura para extrapolar para qsi = 1 e qsi = -1
    Ts[el,:] = [dot(Ne,Tn) dot(Nd,Tn)]; #calcula o valor de T em qsi = 1 e qsi = -1
  end
   Tm = (Ts+[Ts[vs[:,1],2] Ts[vs[:,2],1] ])/2 #calcula o salto entre as temperaturas


  
  for i = 1:ndelinhas_ELEM
    tipo_elem=ELEM[i,end]
    t=range(-1,stop=1,length=tipo_elem);
    td=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem);
    A=zeros(tipo_elem,tipo_elem)
    Al=zeros(tipo_elem-2,tipo_elem)
    A[1,:] = calc_fforma_gen(td[1],tipo,0); #funcao de forma para o no adjacente a esquerda
    A[end,:]= calc_fforma_gen(td[end],tipo,0); #funcao de forma para o no adjacente a direita
    for k=1:tipo_elem-2
      A[k+1,:] = calc_fforma_gen(td[k+1],tipo,0); #funcao de forma para o no adjacente a esquerda
      Al[k,:] = calc_fforma_gen(t[k+1],tipo,dad.afasta); #funcao de forma para o no adjacente a esquerda

    end

    Tsuave[dad.ELEM[i,2:end-1]]=A*[Tm[i,1]; Al*T[dad.ELEM[i,2:end-1]]; Tm[i,2]]
    # @infiltrate
  end
Tsuave
end

