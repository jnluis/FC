%% Problema 2.3 - Oscilador Quártico - Método de Crank-Nicolson
clc;
close all;
clear all;

x0 = 1; % posição inicial
vx0 = 1;
alpha = -0.1;

K  = 1; %N/m
m = 1; % 1 kg

w = sqrt(K/m);
w2 = K/m;

t0= 0;
tf= 20;
h= 0.02;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

vx = zeros(N,1);
vx(1) = vx0;

const = [h/2, K*h/(2*m), 2*alpha];
options = optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

tic
for k=1:N-1
    func = @(xv) fcr(xv,x(k),vx(k),const);
    xv0 = [x(k),vx(k)];
    aux = fsolve(func,xv0,options);
    x(k+1) = aux(1);
    vx(k+1) = aux(2);
end
toc

Em = 0.5*(K*(x.^2+alpha*x.^4)+ m*vx.^2);
subplot(2,2,1)
plot(t,x)
title('Posições ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Posição (m)")

subplot(2,2,2)
plot(t,vx)
title('Velocidades ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Velocidade (m/s)")

subplot(2,2,3)
plot(x,vx)
xlabel("Posição (m)")
ylabel("Velocidade (m/s)")

subplot(2,2,4)
plot(t,Em)
title('Energia mecânica ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Em (J)")

% Localiza máximos
imax = 0;
for k= 2:N-1 % Começa a contar no 2, para ser 0 o indíce 1, acho eu
    if and(x(k+1)-x(k) <=0 , x(k) - x(k-1) >= 0 )
        imax = imax+1;
        aux = lagr(t(k-1:k+1), x(k-1:k+1));
        tmax(imax) = aux(1);
        xmax(imax) = aux(2);
    end
end

nmax= imax;
plf =  polyfit(1:nmax, tmax, 1);
T= plf(1);
A = mean(xmax);

% Solução analítica
Asa = sqrt(x0^2 +m *vx0^2/K);
Tsa = 2*pi/w;

%Mostrar resultados
% Alguns estão comentados porque não são necessários para esta questão
fprintf(' Período método numérico:   %d \n', T);
%fprintf(' Período sol. analítica:    %d \n', Tsa);
%fprintf(' Um a dividir pelo outro:   %d \n', T/Tsa);
fprintf(' Amplitude método numérico: %d \n', A);
%fprintf(' Amplitude sol. analítica:  %d \n', Asa);
%fprintf(' Uma a dividir pelo outra:  %d \n', A/Asa);