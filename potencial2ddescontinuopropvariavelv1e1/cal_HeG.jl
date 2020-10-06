function cal_HeG(dad,npg=8)   # (NOS,ELEM,kmat,afasta,npg=8)
  n = size(dad.NOS,1)         # Quantidade de nos fisicos discretizando o contorno
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  H=zeros(n,n)
  G=zeros(n,n)
  qsi,w = gausslegendre(npg)                     # Quadratura de gauss
  for j=1:nelem                                  # Laço dos elementos (todos os elementos)
    tipo_elem = dad.ELEM[j,end]                  # Tipo do elemento fisico
    x = dad.NOS[dad.ELEM[j,2:1+tipo_elem],:]     # Coordenada (x,y) dos nós geométricos

    for i=1:n                                    # Laço dos pontos fontes (todos os nós)
      pf=dad.NOS[i,:]                            # coordenada pf 
      h = zeros(tipo_elem)                       # array do tamanho "tipo_elem"
      g = zeros(tipo_elem)                       
      nosing = dad.ELEM[j,2:1+tipo_elem] .== i   # Verifica se algum dos nos é ponto fonte
      if sum(nosing) == 1                        # IF elemento contem o ponto fonte
        qsis=range(-1+dad.afasta,stop=1-dad.afasta,length=tipo_elem) # Parametrização de -1 a 1
        eet=qsis[nosing]   #Determina a posição do pf no elemento    # do elemento com a distribuição dos nos
        eet=eet[1]     # transformando em float
        eta,Jt=telles(qsi,eet);
      else       # elemento Ñ contem o pf
        # rel=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        # eet=2*dot(rel,pf-x[1,:])/norm(rel)^2-1     # Ñ ENTENDI, ver: telles
        eets=zeros(tipo_elem-1)
        bs=zeros(tipo_elem-1)
        for ii=1:tipo_elem-1                      # IF elem não contem ponto forte
          rel=x[ii+1,:]-x[ii,:]
          eets[ii]=2*dot(rel,pf-x[ii,:])/norm(rel)^2-1
          N=calc_fforma_gen(eets[ii],tipo_elem,dad.afasta)       #pq calcula isso aqui se embaixo calcula novamente?
          ps=N*x             
          bs[ii]=norm(ps'-pf)                                    # não entendi
        end
        b,ind=findmin(bs)
        eet=eets[ind]                                            # não entendi
        eta,Jt=sinhtrans(qsi,eet,b)
      end
      # qsi,w = otimiza_npg(x,pf)
      # eta,Jt= telles(qsi,eet);
      for k = 1:npg                                   #for para cada ponto de gauss no elemento for cada pf
        # N_geo= calc_fforma_gen(eta[k],tipo_elem_geo,0)#funções de forma (geometrico)
        # N= calc_fforma_gen(eta[k],tipo_elem,afasta)
        # N_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0)#funções de forma generalizada
        N = calc_fforma_gen(eta[k],tipo_elem,dad.afasta)      # N1,N2 e N3
        pg = N*x                                              # Ponto de gauss interpolador
        r = pg'-pf                                            # r= Distancia entre

        # dN_geo = calc_fforma_gen(eta[k],tipo_elem_geo,0)# calcula dN\dξ N1,N2 e N3
        dN = calc_dfforma_gen(eta[k],tipo_elem,dad.afasta)            # calcula dN\dξ 
        dxdqsi = dN*x                                                 # dx/dξ & dy/dξ
        dgamadqsi =norm(dxdqsi)                                       # dΓ/dξ = J(ξ) Jacobiano

        sx=dxdqsi[1]/dgamadqsi;                                       # vetor tangente dx/dΓ
        sy=dxdqsi[2]/dgamadqsi;                                       # vetor tangente dy/dΓ

        kponto=kf(pf)                                                 #definir kno ponto fonte

        Qast,Tast=calsolfund(r,kponto,[sy,-sx])                        # T e Q fundamentais
        # sum(isnan.(Qast)) >= 1 ? println([i,j,k]) : 1 

        h+=N'*Qast*dgamadqsi*w[k]*Jt[k]*kponto
        g+=N'*Tast*dgamadqsi*w[k]*Jt[k]*kponto
      end
      # @show h
      H[i,dad.ELEM[j,2:1+tipo_elem]] += h
      G[i,dad.ELEM[j,2:1+tipo_elem]] += g
    end
  end
  #for i=1:n #Laço dos pontos fontes
    # H[i,i]=0
    # H[i,i]=-sum(H[i,:]) # hipotese de corpo com temperatura constante
    #H[i,i]=-.5
  #end
  H,G
end

function calsolfund(r,kmat,n)
  R=sqrt.(sum(r.^2))
  Qast=sum(r.*n)/R^2/(2*π)       # Equação 4.36
  Tast=-log(R)/(2*π*kmat)
  Qast,Tast
end

function calsolfund_dx_dy(r, kmat, n)
    R = norm(r)
    dQastdx = (n[1]*(r[1]^2 - r[2]^2) + 2*n[2]*r[1]*r[2])/(2*π*R^4)
    dTastdx = r[1]/R^2/(2*π*kmat)
    dQastdy = (n[2]*(-r[1]^2 + r[2]^2) + 2*n[1]*r[1]*r[2])/(2*π*R^4)
    dTastdy = r[2]/R^2/(2*π*kmat)
    dQastdx, dTastdx, dQastdy, dTastdy      # Return
end

function calsolfund_d2x_d2y(r, kmat, n)
    R = norm(r)
    d2Qastdx = (2*r[1]*(n[1]*(r[1]^2 - r[2]^2) + 2*n[2]*r[1]*r[2]) - (n[1]*r[1] + n[2]*r[2])*R^2)/(π*R^6)
    d2Tastdx = (2*r[1]^2 - R^2)/(2*π*kmat*R^4)
    d2Qastdy = (2*r[2]*(n[2]*(r[2]^2 - r[1]^2) + 2*n[1]*r[1]*r[2]) - (n[1]*r[1] + n[2]*r[2])*R^2)/(π*R^6)
    d2Tastdy = (2*r[2]^2 - R^2)/(2*π*kmat*R^4)
    d2Qastdx, d2Tastdx, d2Qastdy, d2Tastdy      # Return
end

function calsolfund_dxdy(r, kmat, n)
  R = norm(r)
  dQdxdy =  -(n[1]*r[2]^3-3*n[2]*r[1]*r[2]^2-3*n[1]*r[1]^2*r[2]+n[2]*r[1]^3)/(π*R^6)
  dTdxdy = r[1]*r[2]/(π*kmat*R^4)
  dQdxdy,dTdxdy
end

function geraCDC(CCSeg,ELEM,MALHA)
  nelem = size(ELEM,1)
  nseg  = size(CCSeg,1)
  tipo_elem=ELEM[1,end]
  tipoCDC=zeros(nelem*tipo_elem)
  valorCDC=zeros(nelem*tipo_elem)
  let
    k=1
    for i=1:nseg  # laço dos segmentos
      for j=1:MALHA[i,end] # laço dos elementos
        for j1=1:tipo_elem # laço dos nós internos ao elemento
         tipoCDC[k] =CCSeg[i,2]
         valorCDC[k] =CCSeg[i,3]
         k+=1
       end
      end
    end
  end
  tipoCDC,valorCDC
end

function aplicaCDC(H,G,dad)
  n=length(dad.tipoCDC)
  A=1*H
  B=1*G
  #temperatura conhecida
  #@infiltrate
  (A[:,1:n])[:,dad.tipoCDC.==0]=  -G[:,dad.tipoCDC.==0]
  B[:,dad.tipoCDC.==0]=  -(H[:,1:n])[:,dad.tipoCDC.==0]

  #fluxo conhecido
  #(A[:,1:n])[:,dad.tipoCDC.==1]=  H[:,dad.tipoCDC.==1]
  #B[:,dad.tipoCDC.==1]=  (G[:,1:n])[:,dad.tipoCDC.==1]

  b=B*dad.valorCDC
  A,b
end

"function aplicaCDC(H,G,dad)
  n=length(dad.tipoCDC)
  A=0*H
  B=0*H
  #temperatura conhecida
  # @infiltrate
  A[:,dad.tipoCDC.==0]=  -G[:,dad.tipoCDC.==0]
  B[:,dad.tipoCDC.==0]=  -H[:,dad.tipoCDC.==0]

  #fluxo conhecido
  A[:,dad.tipoCDC.==1]=  H[:,dad.tipoCDC.==1]
  B[:,dad.tipoCDC.==1]=  G[:,dad.tipoCDC.==1]

  b=B*dad.valorCDC
  A,b
end"

function corrigeCDC(tipoCDC,valorCDC,ELEM,MALHA,NOS,prob)
if prob == "dad_1"
  nelem = size(ELEM,1)
  nseg  = size(CCSeg,1)
  tipo_elem=ELEM[1,end]
  Ta = zeros(nelem*tipo_elem)
  X0 = [0 0] #referencia do sistema de coordenada
  let
    k=1
    for i=1:nseg  # laço dos segmentos
      for j=1:MALHA[i,end] # laço dos elementos
        for j1=1:tipo_elem # laço dos nós internos ao elemento
          if i == 2
            # tipoCDC[k] =CCSeg[i,2]
            r = hypot(NOS[k,1],NOS[k,2])
            theta = acos(NOS[k,1]/r)
            valorCDC[k] = -1/(2*sqrt(r))*(cos(theta/2)*cos(theta) + sin(theta/2)*sin(theta))
          end
          if i == 3
            # tipoCDC[k] =CCSeg[i,2]
            r = hypot(NOS[k,1],NOS[k,2])
            theta = acos(NOS[k,1]/r)
            valorCDC[k] = -1/(2*sqrt(r))*(cos(theta/2)*sin(theta) - sin(theta/2)*cos(theta))
          end
          if i == 4
            # tipoCDC[k] =CCSeg[i,2]
            r = hypot(NOS[k,1],NOS[k,2])
            theta = acos(NOS[k,1]/r)
            valorCDC[k] = 1/(2*sqrt(r))*(cos(theta/2)*cos(theta) + sin(theta/2)*sin(theta))
          end
          r = hypot(NOS[k,1],NOS[k,2])
          theta = acos(NOS[k,1]/r)
          Ta[k] = sqrt(r)*cos(theta/2)
          k+=1
       end
      end
    end
  end
  tipoCDC,valorCDC,Ta
end
end

# function otimiza_npg(x,pf)
#   L = norm(x[end,:]-x[1,:])   # distancia entre extremidades do elemento
#   # Lreal = integra(dgamadqsi,w[k],npg) # comprimento real do elemento
#   r = x .-pf'                 # distancia entre o ponto fonte e todos os nós GEOMÉTRICOS
#   rmin = minimum(sqrt.(r[:,1].^2 + r[:,2].^2))  # obtendo a distancia mínima sqrt(x^2 + y^2) :: Float
#   z = rmin/L
#   if z <= 1
#     # npg = round(Int,K/z)
#   else
#     npg = 4
#   end
#     qsi,w=gausslegendre(npg)
#     qsi1=qsi/2*(a+1).+(a-1)/2
#     qsi2=qsi/2*(1-a).+(a+1)/2
#     qsis=[qsi1; qsi2]
#     ws=[w*(a+1)/2; w/2*(1-a)]
#     qsis,ws
# end
