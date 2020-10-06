## Início da análise
# Para definição do problema escreva em dad, caso a condição de contorno seja dependente
# da posição, ajuste a função corrigeCDC para

# include("potencial.jl")

include("includes.jl")
nelem = 5  #Numero de elementos
npint = 2 #pontos internos
tipo=3       #2:linear ; 3:quadratico; 4:cúbico  (DESCONTINUO)
afasta=1/tipo
npg = 16    #apenas números pares

# PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_0(nelem,npint)
PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_1(nelem,npint)
# PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_2(nelem,npint)   #DT problem
# PONTOS,SEGMENTOS,MALHA,CCSeg,NPX,NPY,k = dad_placa_com_furos(4,20,4)


## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
# afasta=0.5
NOS,ELEM = format_dad(PONTOS,SEGMENTOS,MALHA,tipo,afasta)       # Físicos
NOS,ELEM,vs = format_dad(PONTOS,SEGMENTOS,MALHA,tipo,afasta,true)       # Físicos
tipoCDC,valorCDC = geraCDC(CCSeg,ELEM,MALHA)
dad=dados(NOS,ELEM,tipoCDC,valorCDC,afasta,k,kf) # dado para organizar os argumentos
tipoCDC,valorCDC,Ta = corrigeCDC(tipoCDC,valorCDC,ELEM,MALHA,NOS,"dad_1")   #Ta: temperatura analítica

println("2. Montando a matriz A e o vetor b")
H,G = cal_HeG(dad,npg)  #importante
A,b = aplicaCDC(H,G,dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A\b
println("4. Separando fluxo e temperatura")
T,q = monta_Teq(dad,x) #importante
# _____________________________________________________________________
Tsuave=suavizar(dad,vs,T)
## Calculo das variaveis no contorno e em pontos internos ______________
# println("5. Gerando pontos internos")
# PONTOS_INT = gera_p_in(NPX,NPY,PONTOS,SEGMENTOS); # gera os pontos internos
# # ______________________________________________________________________

# println("6. Calculando temperatura e fluxo nos pontos internos")
# Ti = domain_T(PONTOS_INT,dad,T,q,npg) #(...,npg) Calcula a T e q nos Pt internos
# dTidx, dTidy,d2Tidx,d2Tidy,dTidxdy=d_SF(PONTOS_INT,dad,T,q,8)
# Ti_ANA = dad1(PONTOS_INT)

# Q_ANA = Qin_analitic(PONTOS_INT) #dad_1
# Qc_ANA = Qin_analitic(NOS) #dad_1

# ## Calcula a  derivada da temperatura no contorno _________
# println("7. Calculando derivada da temperatura no contorno")

# dT = cal_dTSF(NOS,dad,T,q)

# dT2 = cal_d2TSF(NOS,dad,T,q)