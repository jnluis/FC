%% Problema 4.3 - Método de diferenças finitas — perfil de temperaturas numa resistência elétrica cilíndrica
clc;
clear all;
close all;

ordem2=false;
ordem2=true;

% Para usar a aproximação de 1ºordem para a cond. fronteira de Neumann
% comenta-se a linha 7
% Para usar a aproximação de 2ºordem para a cond. fronteira de Neumann
% comenta-se a linha 6

Tamb=20;
lambda=0.1;
Q=2.1e6;
r0=0;
R=0.001;

h=1e-6;
r=r0:h:R;
N=numel(r);

T=zeros(N,1);
A=-2*eye(N);
% cria uma matriz identidade de tamanho N*N
% e multiplica os elementos por -2

A(N,N)=1;

for i=2:N-1
    A(i,i+1)=1+h/(2*r(i)); % diagonal acima da principal
    A(i,i-1)=1-h/(2*r(i)); % diagonal acima da principal
end

b=(-h^2*Q/lambda)*ones(N,1); % isto é o que está do lado direito da equação
b(N)=Tamb;

if ordem2
    A(1,1)=-2;
    A(1,2)=2;
else
    A(1,1)=-1;
    A(1,2)=1;
    b(1)=0;
end

T= linsolve(A,b);
plot(r,T);
xlabel('{\it r}');
ylabel('{\it T}');
[maxT, indice]=max(T);

fprintf('\nTmax = %f, em r =%f\n\n',maxT,r(indice))