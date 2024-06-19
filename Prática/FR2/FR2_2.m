%%% Problema FR2.2: BVP – Método do Shooting
% Feito com Euler-Cromer em vez de ode45
clc;
clear all;
close all;

h= 0.001;
tf=50;
t = 0:h:tf;
N=length(t);

km=2;
m=1.5;
alpha=-0.2;
x0=1.9;
v0=0;

x=zeros(N,1);
v=zeros(N,1);

x(1)=x0;
v(1)=v0;

for k=1:N-1
    % Euler-Cromer
    v(k+1)=v(k)+( (-km/m) * x(k)*(1+(3/2)*alpha*x(k)) ) *h;
    x(k+1)=x(k)+v(k+1)*h;
end

figure(1)
plot(t,x)
xlabel('t');
ylabel('x');