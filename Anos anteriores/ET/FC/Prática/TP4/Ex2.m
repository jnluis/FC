%%  Ex4.2

clear,clc,close all

%% CONSTANTES

u = 10e-3;
L= 1;
T = 10e3;
n=1;
h=0.001;

x=0:h:L;
N=length(x);

%solução analitica da frequencia
w1= n*pi/L*sqrt(T/u);

%% METODO DAS DIFERENÇAS FINITAS

A = eye(N-2,N-2);   %criar uma matriz identidade com N-2 linhas e colunas
A = -2*A;           %diagonal = -2

for k=1:N-2
    A(k,k+1)= 1;    %fazer o triangulo superior
    A(k+1,k)= 1;    %fazer o triangulo inferior
end    

%% CALCULO DOS VALORES PROPRIOS

[vec,val]=eigs(A,3,'sm')    %calculo dos valores proprios

w=diag( sqrt(-(val*T)/(u*h^2)) ) %calculo dos omegas para os diferentes 
                                %fazer expressão em função da frequencia

