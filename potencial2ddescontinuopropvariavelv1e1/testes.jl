# include("potencial.jl")
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

            tipo_elem_geo = dad.ELEM_GEO[j,end]   #Tipo do elemento geometrico
            tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
            hx,hy = zeros(tipo_elem),zeros(tipo_elem)
            gx,gy = zeros(tipo_elem),zeros(tipo_elem)    # Array de 'tipo_elem' linhas
      
            x = dad.NOS_GEO[dad.ELEM_GEO[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
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
                N_geo=calc_fforma(eet,tipo_elem_geo,0)
                ps=N_geo*x
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

            tipo_elem_geo = dad.ELEM_GEO[j,end]   #Tipo do elemento geometrico
            tipo_elem = dad.ELEM[j,end]           #Tipo do elemento fisico
            hxx,hyy,hxy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)
            gxx,gyy,gxy = zeros(tipo_elem),zeros(tipo_elem),zeros(tipo_elem)    # Array de 'tipo_elem' linhas
      
            x = dad.NOS_GEO[dad.ELEM_GEO[j,2:1+tipo_elem],:]   # Coordenada (x,y) dos nós geométricos
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
                N_geo=calc_fforma(eet,tipo_elem_geo,0)
                ps=N_geo*x
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
# dT = cal_dTSF(NOS,dad,T,q)

# dT2 = cal_d2TSF(NOS,dad,T,q)
# dTa=Qin_analitic(NOS);
# T_ANA = dad1(NOS);

# qa=dTa[:,1];qa[1:3*nelem]=-dTa[1:3*nelem,2];qa[6*nelem+1:9*nelem]=dTa[6*nelem+1:9*nelem,2];qa[15*nelem+1:end]=-dTa[15*nelem+1:end,2];
# qa[9*nelem+1:12*nelem]=dTa[9*nelem+1:12*nelem,2];qa[12*nelem+1:15*nelem]=-qa[12*nelem+1:15*nelem];
# dT1 = cal_dTSF(NOS,dad,T_ANA,qa);

# [q dTa dT]
# plot(dT1+dTa)
# plot(hxv)
# plot(hyv)
# plot(gxv)
# plot(gyv)
