function derivadas=f_FR2_2(t,solucao,K,m,alfa)

derivadas=zeros(2,1);
x=solucao(1);
v=solucao(2);

derivadas(1)=v;
derivadas(2)=-K/m * x *(1+3/2*alfa*x);