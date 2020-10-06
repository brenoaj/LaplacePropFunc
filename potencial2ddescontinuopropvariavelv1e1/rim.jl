function Monta_M_RIM(NOS, ELEM, PONTOS_INT, k, npg1, npg2)
    n_nos = size(NOS, 1)
    nelem = size(ELEM, 1)
    n_noi = length(PONTOS_INT[:,1]); #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi;
    nodes = [NOS;PONTOS_INT]
    M = zeros(n_pontos, n_pontos);
    F = zeros(n_pontos, n_pontos);

# Cálculo da matriz [F]
    for i = 1:n_pontos
        xi = nodes[i,1];
        yi = nodes[i,2];
        for j = 1:n_pontos
            xj = nodes[j,1];
            yj = nodes[j,2];
            r = sqrt((xi - xj)^2 + (yi - yj)^2);
            F[i,j] = interpola(r);       
        end
    end
    for i = 1:n_pontos #Laço dos pontos radiais
        pr = nodes[i,:]
        for j = 1:n_pontos #Laço dos pontos fontes
            pf = nodes[j,:]
            for el = 1:nelem
                tipo = ELEM[el,5]; # Tipo de elemento (1=cont�nuo e 2=descont�nuo)
                x1 = nodes[ELEM[el,2],:]
                x2 = nodes[ELEM[el,3],:]
                x3 = nodes[ELEM[el,4],:]
                 m_el = calc_m(x1, x2, x3, pr, pf, k, qsi1, w1, npg2, tipo)
                M[j,i] = M[j,i] + m_el
            end
        end
    end
    M / F
end

function interpola(r)
    	if r == 0
        return 0
    end
    r^2 * log(r)
    # r

end
function interpola(r, i)
    (4 * r^4 * log(r) - r^4) / 16
    # r^3/3
end
function  calc_m(x1, x2, x3, pr, pf, k, qsi1, w1, npg2, tipo)
    npg = length(w1);
    m_el = 0;

    for i = 1:npg
        N1, N2, N3 = calc_fforma(qsi1[i], tipo); # Calcula as fun��es de forma
        dN1dqsi, dN2dqsi, dN3dqsi = calc_dfforma(qsi1[i], tipo); # Calcula as
        # derivadas das fun��es de forma

         pg = N1 * x1 + N2 * x2 + N3 * x3; # Calcula a coordenada x do ponto de integra��o

        dxdqsi = dN1dqsi * x1 + dN2dqsi * x2 + dN3dqsi * x3;
        dgamadqsi = norm(dxdqsi)

        sx = dxdqsi[1] / dgamadqsi; # Component x do vetor tangente
        sy = dxdqsi[2] / dgamadqsi; # Componente y do vetor tangente

        nx = sy; # Componente x do vetor normal unit�rio
        ny = -sx; # Componente y do vetor normal unit�rio

        qsi2, w2 = gausslegendre(npg2)
        r = pg - pf
        m = calcula_F(pr, pf, pg, [nx,ny], k, qsi2, w2);
        
        m_el += dot([nx,ny], r) / norm(r)^2 * m * dgamadqsi * w1[i]
    end
    return m_el
end

function calcula_F(pr, pf, pg, n, k, qsi, w);
    npg = length(w);
    R = (pg - pf)
    r = norm(R)
    drodqsi = r / 2; # Jacobiano da transforma��o de vari�vel de r
           #    para qsi (dr/dqsi)
    F = 0; # Inicializa a integral de F_area
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1); # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        Tast = -log(ro) / (2 * π * k)
        f = interpola(rline)
         F = F + Tast * f * norm(ro) * drodqsi * w[i]; # Integral de F_area
    end
    return F
end


function Monta_M_RIMd(dad, PONTOS_INT, npg1, npg2)          #nome de variáveis modificados 
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(PONTOS_INT,1); #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi;
    nodes = [dad.NOS;PONTOS_INT]
    S = zeros(n_pontos);
    S1 = zeros(n_pontos);
    F = zeros(n_pontos, n_pontos);
    D = zeros(n_pontos, n_pontos);

# Cálculo da matriz [F] e [D]
    for i = 1:n_pontos                #laço dos pontos fonte  
        xi = nodes[i,1];
        yi = nodes[i,2];       
        for j = 1:n_pontos            #laço dos pontos de interpolação                       
            xj = nodes[j,1];
            yj = nodes[j,2];
            r = sqrt((xi - xj)^2 + (yi - yj)^2);
            
            F[i,j] = interpola(r);                  # fm função de interpolção   Matriz F PG Áquila
           
            udx = (xi-xj)/(r*π);                           #d(log(r)/2pi)/dx  = d(r)dx  *  1/(r*2*pi)    era (xi-xj)*π*2/(r^2)
            udy = (yi-yj)/(r*π);                           #d(log(r)/2pi)/dy  = d(r)dy  *  1/(r*2*pi)
            kfdxdy= ForwardDiff.gradient(kf,[xi yi]);      #valores da derivada de kf para x e para y
            D[i,j] = udx*kfdxdy[1] + udy*kfdxdy[2]             #NOVO D : u*(x,qsi),i  *  k(x),i     (eq 8 artigo Loeffler)
            
            #D[i,j] = -log(r) / (2 * π * k)          #u*(qsi,X)  É o resto da eq 4.33 (sem M(s) e F)
        end
    end
    for i = 1:n_pontos #Laço dos pontos radiais
        pr = nodes[i,:]
        #for j = 1:n_pontos #Laço dos pontos fontes
        j=i    
        pf = nodes[j,:]
            
            for el = 1:nelem
                tipo = ELEM[el,end]; # Tipo de elemento (1=cont�nuo e 2=descont�nuo)
                X = nodes[dad.ELEM[el,2:1+tipo],:] 
                s_el, s_el1 = calc_s(X, pr, pf, dad.kf, qsi1, w1, npg2, tipo)
                S[i] = S[i] + s_el                    #calculo da matriz s no PG Áquila
                #S1[j] = S1[j] + s_el1
             end
        #end 
    end
	# @show size(M
	# @show length(M)
    M = ones(length(S)) * (S' / F) .* D          # M é matriz I1
    for i = 1:n_pontos #Laço dos pontos radiais
        M[i,i] = 0
        M[i,i] = -sum(M[i,:])                    # diagonal de I1 Sii
    end
    #M+diagm(0 => S1)                             # I2 diagonal é adicionado a I1
    return M
end

function  calc_s(X, pr, pf, kf, qsi1, w1, npg2, tipo)         #nome modificado
    npg = length(w1);
    s_el, s_el1 = 0, 0

    for i = 1:npg
        N = calc_fforma_gen(qsi1[i], tipo,dad.afasta); # Calcula as fun��es de forma
        dNdqsi = calc_dfforma_gen(qsi1[i], tipo,dad.afasta); # Calcula as derivadas das fun��es de forma
        pg = N*X    # Ponto de gauss interpolador

        dXdqsi = dNdqsi*X   # [dx/dξ dy/dξ]
        dgamadqsi = norm(dXdqsi)   #Jacobiano

        sx = dXdqsi[1] / dgamadqsi; # Component x do vetor tangente
        sy = dXdqsi[2] / dgamadqsi; # Componente y do vetor tangente

        nx = sy; # Componente x do vetor normal unit�rio
        ny = -sx; # Componente y do vetor normal unit�rio
        qsi2, w2 = gausslegendre(npg2)
        # @infiltrate
        r, r1 = pg' - pr, pg' - pf
        R, R1 = norm(r), norm(r1)
        s = interpola(R, 1)
        kponto=kf(pf)  
        #s1 = -(2 * R1^2 * log(R1) - R1^2) / 4 / (2 * π * kponto)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        s_el += dot([nx,ny], r) * s * dgamadqsi * w1[i] / (norm(r)^2)        #eq. s=... apos 4.33  Áquila   norm(r)^2 n deve ser ^2
        #s_el1 += dot([nx,ny], r1) / norm(r1)^2 * s1 * dgamadqsi * w1[i]
    end
    return s_el, s_el1
end

function calcula_Fd(pr, pf, pg, n, k, qsi, w);

    npg = length(w);
    R = (pg - pf)
    r = norm(R)
    drodqsi = r / 2; # Jacobiano da transforma��o de vari�vel de r
         #    para qsi (dr/dqsi)
    F, F1 = (0, 0); # Inicializa a integral de F_area
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1); # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        Qast, Tast = calsolfund(xc - pf, k, n)
        f = interpola(rline)
        F = F + f * norm(ro) * drodqsi * w[i]; # Integral de F_area
        F1 = F1 + Tast * norm(ro) * drodqsi * w[i]; # Integral de F_area
    end
    return F, F1
end