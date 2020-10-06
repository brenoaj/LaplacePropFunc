function calc_fforma(qsi,tipo,a)
if tipo==2 # Elemento linear
  N1 = (-qsi+(1-a))./(2-2*a)  # Função de forma N1
  N2 = (qsi+(1-a))./(2-2*a)   # Função de forma N2
  N= [N1 N2]
elseif tipo==3  # Elemento quadratico
    N1 = qsi.*(qsi-(1-a))/(-(1-a)-(1-a))/(-(1-a)); # Função de forma N1 => quadrática descontínua
    N2 = (qsi+(1-a)).*(qsi-(1-a))/((1-a))/(-(1-a)); # Função de forma N2 => quadrática descontínua
    N3 = qsi.*(qsi+(1-a))/((1-a)+(1-a))/((1-a));  # Função de forma N3 => quadrática descontínua
    N= [N1 N2 N3]
end
N
end

function calc_fforma_gen(qsi,tipo,a) # funções de forma gerais para m elementos
  dN = zeros(1,tipo) #vetor linha das funções de forma N_1, N_2, ..., N_m
  N = ones(typeof(qsi),1,tipo) #vetor linha das funções de forma N_1, N_2, ..., N_m
  qsis=range(-1+a,stop=1-a,length=tipo) # Parametrização de -1 a 1 ξ(x_i)
  for i = 1:tipo
    for j = 1:tipo  # m(numero de nós) = tipo
      if j != i
        N[i] = N[i]*(qsi - qsis[j])/(qsis[i] - qsis[j]) #Eq. interpolador de lagrange
      end
    end
  end
  N
end

function calc_dfforma_gen(qsi,tipo,a) # funções de forma gerais para m elementos
  ForwardDiff.derivative(x->calc_fforma_gen(x,tipo,a),qsi)
end

function calc_dfforma(qsi,tipo,a)
  if tipo==2 # Elemento linear
    N1 = -1.0/(2-2*a)# forma N1
    N2 = 1.0/(2-2*a)# forma N3
    N= [N1 N2]
  elseif tipo==3
   # Elemento quadratico
      N1 = (qsi-(1-a))/(-(1-a)-(1-a))/(-(1-a))+(qsi)/(-(1-a)-(1-a))/(-(1-a));# Função de forma N1 => quadrática descontínua
      N2 = (qsi-(1-a))/((1-a))/(-(1-a))+(qsi+(1-a))/((1-a))/(-(1-a))# Função de forma N2 => quadrática descontínua
      N3 = (qsi+(1-a))/((1-a)+(1-a))/((1-a))+qsi/((1-a)+(1-a))/((1-a));# Função de forma N3 => quadrática descontínua
      N= [N1 N2 N3]
  end
  N
end

function calc_dfforma1(qsi,tipo)
  if tipo==1
    dN1dqsi=-0.5 + qsi;
    dN2dqsi=-2*qsi;
    dN3dqsi=0.5 + 1.0*qsi;
  elseif tipo==2
    dN1dqsi=3*(-1 + 3*qsi)/4;
    dN2dqsi=-9*qsi/2;
    dN3dqsi=3*(1 + 3*qsi)/4;
  elseif tipo==3
    dN1dqsi=-.5+qsi*0;
    dN2dqsi=0+qsi*0;
    dN3dqsi=.5+qsi*0;
  end
  dN1dqsi,dN2dqsi,dN3dqsi
end
