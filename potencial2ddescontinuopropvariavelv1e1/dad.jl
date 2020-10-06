# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function problema1d(ne=15,npi = 11)
    PONTOS  = [1 0 0
        2 1 0
        3 1 1
        4 0 1 ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [1 ne
        2 ne
        3 ne
        4 ne]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [1 1 0
        2 0 1
        3 1 0
        4 0 0]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    NPX = npi
    NPY = npi
    return PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k
end

function dad_1(ne=15,npi = 11)
    # Matriz para definição de pontos que definem a geometria
    # PONTOS = [número do ponto, coord. x do ponto, coord. y do ponto]
    PONTOS  = [1 0 0
        2 1 0
        3 1 1
        4 0 1]
    # Segmentos que definem a geometria
    # SEGMENTOS=[N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    #= Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
                              inicial para o ponto final)
                         < 0 -> O centro é à direita do segmento (do ponto
                              inicial para o ponto final)
                         = 0 -> O segmento é uma linha reta        =#
    SEGMENTOS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0]
    # Matriz para definição da malha
    # MALHA = [n° do segmento, n° de elementos no segmento]
    MALHA = [1 ne
        2 ne
        3 ne
        4 ne]
    #= Condições de contorno nos segmentos
       CCSeg=[no do segmento, tipo da CDC, valor da CDC]
       tipo da CDC = 0 => a temperatura é conhecida
       tipo da CDC = 1 => O fluxo é conhecido         =#
    CCSeg=[1 1 0
        2 1 -1
        3 1 0
        4 0 0]
    k = 1                                # Condutividade Térmica do material
    kf = f(x) = x[1]*1 + x[2]*0 + 1        # função kmat variável

    NPX = npi
    NPY = npi
    return PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k,kf
end

function dad_2(ne=15,npi = 11)  # Problema para otimização topológica
    # Matriz para definição de pontos que definem a geometria
    # PONTOS = [número do ponto, coord. x do ponto, coord. y do ponto]
    PONTOS = [1  0  0
              2 0.2 0
              3 0.8 0
              4  1  0
              5  1 0.2
              6  1  1
              7 0.6 1
              8 0.4 1
              9  0  1
              10 0 0.2]
    # Segmentos que definem a geometria
    # SEGMENTOS=[N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    #= Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
                              inicial para o ponto final)
                         < 0 -> O centro é à direita do segmento (do ponto
                              inicial para o ponto final)
                         = 0 -> O segmento é uma linha reta        =#
    SEGMENTOS= [1 1 2 0
                2 2 3 0
                3 3 4 0
                4 4 5 0
                5 5 6 0
                6 6 7 0
                7 7 8 0
                8 8 9 0
                9 9 10 0
               10 10 1 0]
    # Matriz para definição da malha
    # MALHA = [n° do segmento, n° de elementos no segmento]
    MALHA = [1 ne
             2 ne
             3 ne
             4 ne
             5 ne
             6 ne
             7 ne
             8 ne
             9 ne
            10 ne]
    #= Condições de contorno nos segmentos
       CCSeg=[no do segmento, tipo da CDC, valor da CDC]
       tipo da CDC = 0 => a temperatura é conhecida
       tipo da CDC = 1 => O fluxo é conhecido         =#
       Th = 350
       Tl = 300
    CCSeg = [1 0 Th
             2 1 0
             3 0 Th
             4 0 Th
             5 1 0
             6 1 0
             7 0 Tl
             8 1 0
             9 1 0
            10 0 Th]
    k = 1     # Condutividade Térmica do material
    NPX = npi
    NPY = npi
    return PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k
end

function dad_placa_com_furos(nfuros,ne_aresta,ne_furos)
L=1; # Lado da placa (placa quadrada)
Area=12.47*L^2/(100*nfuros^2); # Área dos furos = 12.47# da área da placa
raio=√(Area/π)  ; # Raio dos furos (diminui com o aumento dos furos)
deltax=L/nfuros; # Distância entre os centros dos furos na direção x
deltay=deltax; # Distância entre os centros dos furos na direção y
cy=deltay/2;
ifuro=1; # Inicializa o contador de furos
k=1; # Condutividade térmica a placa

# NPX=13; # Números de pontos internos na direção x
# NPY=13; # Número de pontos internos na direção y
NPX=11; # Números de pontos internos na direção x
NPY=11; # Número de pontos internos na direção y


PONTOS=zeros(4*nfuros^2+4,3); # Inicializa a matriz PONTOS
SEGMENTOS=zeros(4*nfuros^2+4,4); # Inicializa a matriz PONTOS

MALHA=zeros(UInt32,4*nfuros^2+4,2); # Inicializa a matriz MALHA

CCSeg=zeros(4*nfuros^2+4,3);# Inicializa a matriz CCSeg

# Coordenada dos pontos que definem o contorno externo da placa
# PONTO = [número do ponto, coord. x do ponto, coord. y do ponto]
PONTOS[1:4,1:3]=[1 0 0; 2 L 0;3 L L;4 0 L]

# Segmentos que definem o contorno externo da placa
#  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
#                                                  Radio, tipo do elemento]
# Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
#                          inicial para o ponto final)
#                   < 0 -> O centro é à direita do segmento (do ponto
#                          inicial para o ponto final)
#                   = 0 -> O segmento é uma linha reta
SEGMENTOS[1:4,1:4]=[1 1 2 0
    2 2 3 0
    3 3 4 0
    4 4 1 0];

# Quantidade de elementos no contorno externo da placa
# MALHA = [no do segmento, no de elementos no segmento]
MALHA[1:4,:]=[1 ne_aresta
    2 ne_aresta
    3 ne_aresta
    4 ne_aresta];

# Condições de contorno nos segmentos
# CCSeg=[no do segmento, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => a temperatura é conhecida
# tipo da CDC = 1 => O fluxo é conhecido
CCSeg[1:4,:]=[1 1 0
    2 0 1
    3 1 0
    4 0 0];

# Cria os furos e aplica as condições de contorno
for i=1:nfuros
    cx=deltax/2;
    for j=1:nfuros
        # Coordenada dos pontos que definem os furos
        # Os pontos que definem os furos são dois pontos diametralmente
        # opostos, com a coordenada y = constante.

        PONTOS[4*ifuro+1:4*ifuro+4,:]=[4*ifuro+1 cx-raio cy;
        4*ifuro+2 cx cy+raio;
        4*ifuro+3 cx+raio cy;
        4*ifuro+4 cx cy-raio]

        # Segmentos que definem os furos
        # Os segementos são definidos pelos pontos diametralmente opostos e
        # pelo raio dos furos
        SEGMENTOS[4*ifuro+1:4*ifuro+4,1:4]=[4*ifuro+1 4*ifuro+1 4*ifuro+2 -raio;
            4*ifuro+2 4*ifuro+2 4*ifuro+3 -raio;
            4*ifuro+3 4*ifuro+3 4*ifuro+4 -raio;
            4*ifuro+4 4*ifuro+4 4*ifuro+1 -raio];

        CCSeg[4*ifuro+1:4*ifuro+4,1:3]=[4*ifuro+1 1 0
            4*ifuro+2 1 0
            4*ifuro+3 1 0
            4*ifuro+4 1 0];


        MALHA[4*ifuro+1:4*ifuro+4,1:2]=[4*ifuro-3 ne_furos
            4*ifuro+2 ne_furos
            4*ifuro+3 ne_furos
            4*ifuro+4 ne_furos]

        ifuro=ifuro+1
        cx=cx+deltax
    end
    cy=cy+deltay
end
return PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k
end
