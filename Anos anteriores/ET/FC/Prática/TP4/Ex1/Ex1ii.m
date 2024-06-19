%%  Ex4.1 ii.

clear,clc,close all

%% CONSTANTES

u = 10e-3;
L= 1;
T = 10e3;
n=1;
h=0.001;

x=0:h:L;
N=length(x);

Dy=zeros(1,N);
y=zeros(1,N);

Dy(1)=2*10^-2;
y(1)=0;

%% METODO DE EULER-CROMER

w = 2000;

for k=1:N-1     %APLICAÇÃO DO EULER-CROMER
    Dy(k+1)=Dy(k)-w^2*u/T*y(k)*h;
    y(k+1)=y(k) + Dy(k+1)*h;
end

%% GRÁFICOS

plot(x,y)
title('y em função de x')
xlabel('x')
ylabel('y')

%Neste caso, a frequência utilizada não é uma frequencia fundamental pois
%inicia em 0, mas não acaba em 0

%% CALCULO DE Y(L)

yL=y(L);   %y(L) = 0