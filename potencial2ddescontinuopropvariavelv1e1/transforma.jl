function telles(gamm,eet)

eest = eet^2 - 1;
term1 = eet*eest + abs(eest);
if term1 < 0
    term1 = (-term1)^(1/3);
    term1 = -term1;
else
    term1 = term1^(1/3);
end

term2 = eet*eest - abs(eest);
if term2 < 0
    term2 = (-term2)^(1/3);
    term2 = -term2;
else
    term2 = term2^(1/3);
end
GAMM = term1 + term2 + eet;


Q = 1 + 3*GAMM^2;
A = 1/Q;
B = -3*GAMM/Q;
C = 3*GAMM^2/Q;
D = -B;

eta = A*gamm.^3 + B*gamm.^2 .+ C*gamm .+ D;
Jt = 3*A*gamm.^2 .+ 2*B*gamm .+ C;
return  eta,Jt
end


function sinhtrans(u,a,b)
    #https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.1208
    # https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.2244 iterado
    # https://www.sciencedirect.com/science/article/pii/S0955799715001125 3d
    if b == 0
        b = 1e-6
    end
    μ=1/2*(asinh((1+a)/b)+asinh((1-a)/b))
    η=1/2*(asinh((1+a)/b)-asinh((1-a)/b))

    # x=(asinh.((u.-a)./b).+η)./μ
    x=a.+b*sinh.(μ*u.-η)
    J=b*μ*cosh.(μ*u.-η)
    x,J
end

function Ma(u,a,b)
# https://link.springer.com/content/pdf/10.1007%2Fs00466-002-0340-0.pdf
# recomenda dividir a integral em 2 na singularidade
J=sqrt.(b.^2 .+(u .-a).^2)
# x=log.(J.+(u.-a))
# xi=log.(sqrt.(b.^2 .+(-1 .-a).^2).+(-1. -a))
# xf=log.(sqrt.(b.^2 .+(1 .-a).^2).+(1. - a))
x=1/2*(exp.(u).-b^2*exp.(-u)).+a
xi=1/2*(exp.(-1).-b^2*exp.(1)).+a
xf=1/2*(exp.(1).-b^2*exp.(-1)).+a
x=2*(x.-(xf+xi)/2)/(xf-xi)
x,J
end

function Xie(u,a,b)
# https://www.sciencedirect.com/science/article/pii/S0955799711000208



end
#
# f(x)=(1 .-x.^2)./sqrt.((x.-a).^2 .+b^2)
# f1(x)=.5*log.((x.-a).^2 .+b^2)
# a=1
# b=0.1
# ana=(-0.5 *(a - 1)* (log((a - 1)^2 + b^2) - 2) - b *atan((a - 1)/b))-(-0.5* (a + 1)* (log((a + 1)^2 + b^2) - 2) - b*atan((a + 1)/b))
#
# u,w=gausslegendre(10);
# x,j=sinhtrans(u,a,b+1e-6);
# e,jt=telles(u,a);
# m,jm=Ma(u,a,b);
#
# [sum(f1(u).*w) sum(f1(x).*w.*j) sum(f1(e).*w.*jt) sum(f1(m).*w.*jm)].-ana
