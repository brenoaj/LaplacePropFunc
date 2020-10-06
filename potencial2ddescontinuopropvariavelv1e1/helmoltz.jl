## Início da análise
# Para definição do problema escreva em dad, caso a condição de contorno seja dependente
# da posição, ajuste a função corrigeCDC para

include("includes.jl")

nelem = 5                           # Numero de elementos por segmento
npint = 3                           # pontos internos
tipo= 2                             # 2:linear ; 3:quadratico; 4:cúbico  (DESCONTINUO)
afasta=1/tipo                       # Afastamento do descontínuo
npg = 8                             # NPG: apenas números pares
W=1                                 # frequencia
c=1                                 # velocidade do som
# PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_0(nelem,npint)
PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k,kf = dad_1(nelem,npint)
# PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_2(nelem,npint)   #DT problem
# PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_placa_com_furos(4,20,4)


# afasta=0.5
#NOS,ELEM = format_dad(PONTOS,SEGMENTOS,MALHA,tipo,afasta)       # Físicos
NOS,ELEM,vs = format_dad(PONTOS,SEGMENTOS,MALHA,tipo,afasta,true)       # Físicos
scatter(NOS[:,1],NOS[:,2])
tipoCDC,valorCDC = geraCDC(CCSeg,ELEM,MALHA)
dad=dados(NOS,ELEM,tipoCDC,valorCDC,afasta,k,kf) # dado para organizar os argumentos
#tipoCDC,valorCDC,Ta = corrigeCDC(tipoCDC,valorCDC,ELEM,MALHA,NOS,"dad_1")   #Não é usado nesse problema


H,G = cal_HeG(dad,npg)  #importante
nc=size(H,1)
PONTOS_INT = gera_p_in(NPX,NPY,PONTOS,SEGMENTOS); # gera os pontos internos
Hi,Gi=Calc_HeG_interno(PONTOS_INT,dad,npg,tipo)
ni=size(Hi,1)
Ht=[H zeros(nc,ni);Hi -zeros(ni,ni)] 
Gt=[G;Gi]


include("rim.jl")
M=Monta_M_RIMd(dad,PONTOS_INT,npg,npg)
#Mf=((W/c)^2)*M                             
H2 = Ht + M 
for i = 1:nc                                #i=1:size(dad.NOS,1) #Laço dos pontos fontes
  H2[i,i]=0
  H2[i,i]=-sum(H2[i,:]) # hipotese de corpo com temperatura constante
end


#A1,b1 = aplicaCDC(Gt,Ht+Mf[:,1:nc],dad)
A1,b1 = aplicaCDC(H2,Gt,dad)   # Calcula a matriz A e o vetor b        
x1 = A1\b1                     #é o vetor de condições desconhecidas em cada NÓ 

T,q = monta_Teq(dad,x1)          # T temp e q fluxo do contorno       ok
Tin=x1[nc+1:end]                 # Temp interna                  ok
TinAna= 2 .* log.(1 .+ PONTOS_INT[:,1])
Ein=(Tin.-TinAna)./TinAna
plot(Ein)


TAna = 2 * log.(1 .+ dad.NOS[:,1]) 
plot(T)            
plot!(TAna) 

qAna = 2 ./ ( 1 .+ dad.NOS[31:40])          
#plot(qAna)
plot(q)                                # ok obs: sol fund é para qx, calculamos qxy e CDC são para qx ou qy




#TAin=dad1(PONTOS_INT)            #Temperatura interna analítica  introduzido
#erroTin=(Tin-TAin)/TAin           #erro                          introduzido
 
#Fazer calculos de erros do 1o exemplo do artigo
# Eq de sol analit no Artigo

#QinA=Qin_analitic(PONTOS_INT)    #Fluxo interno analítico        introduzido


#Ti = domain_T(PONTOS_INT,dad,T,q,8,PONTOS_INT)                  #introduzido

#q-QinA                                                          #introduzido