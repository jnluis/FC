%%% Problema FR2.1: BVP – Método do Shooting
% Feito com Euler-Cromer
clc;
clear all;
close all;

h= 0.01;
tf=7;
x = 0:h:tf;
N=length(x);

beta=18;
T=zeros(N,1);
DT=zeros(N,1);

v=1.7187;
DT0=10^-4*v;

T(1)=10^-4;
DT(1)=DT0;

for k=1:N-1
    % Euler-Cromer
    DT(k+1) = DT(k) + (v*DT(k)-beta * exp(-1/T(k))*(1+DT(k)-v*T(k)))*h;
    T(k+1) = T(k)+DT(k+1)*h;

end

figure(1)
subplot(1,2,1)
plot(x,T)
xlabel('x');
ylabel('Temperatura T');

subplot(1,2,2)
plot(x,DT)
xlabel('x');
ylabel('Derivada de T');