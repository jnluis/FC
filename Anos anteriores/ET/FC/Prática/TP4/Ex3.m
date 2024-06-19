%%  Ex3

clear,clc,close all

%%  Constantes

R=1e-3;
T_amb=20;
Q=2.1e6;
lameda=0.1;

h=1e-5;
r=0:h:R;

N=length(r);

%%  Resolução utilizando o linsolve

A=eye(N,N);     %matriz identidade
A=-2*A;         %diagonal igual a -2
A(N,N)=1;       %ultimo valor da diagonal =1

for k=2:N-1
    A(k,k+1)= 1+h/(2*r(k));    %fazer a diagonal superior
    A(k,k-1)= 1-h/(2*r(k));    %fazer a diagonal inferior
end

A(1,2)=2;
b=-h^2*Q/lameda * ones(N,1);    %matriz b
b(N)=T_amb;                     %ultimo valor =20

T=linsolve(A,b);

%%  Gráfico

plot(r,T)
