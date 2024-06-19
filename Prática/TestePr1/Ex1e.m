%% Ex1e
clc;
close all;
clear all;

x0 = 0.1; % em m
vx0 = 0;
y0= 0.1;
vy0 = 0;

K  = 0.4; %N/m
m = 0.2; % 200 g
L = 1;

w = K/m;
w2 = w^2;
g = 9.81;

t0= 0;
tf= 30;
h= 0.001;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

vx = zeros(N,1);
vx(1) = vx0;

y = zeros(N,1);
y(1) = y0;

vy = zeros(N,1);
vy(1) = vy0;

for k=1:N-1
    vx(k+1) = vx(k) + (-g*(x(k)/L) - K*(x(k)-y(k)))*h;
    vy(k+1) = vy(k) + (-g*(y(k)/L) + K*(x(k)-y(k)))*h;
    
    x(k+1) = x(k) + vx(k+1) *h;
    y(k+1) = y(k) + vy(k+1) *h;
end

%         disp("Método de Euler Implícito com linsolve")
%         A = [1 -h; w^2*h 1];
%         for k=1: N-1
%             % Definindo o vetor b com base nas equasções do método de Euler
%             % Implícito
%             b = [x(k); vx(k)];
% 
%             % Resolvendo o sistema linear utilizando linsolve
%             % X = linsolve( A, B) solves the linear system AX = B
%             aux = linsolve(A, b);
% 
%             % Atualizando as variáveis de posição e velocidade
%             vx(k+1) = aux(2);
%             x(k+1) = aux(1);
%         end

figure(1)
plot(t,x)
title('X ao longo do tempo');
xlabel("Tempo (s)")
ylabel("X (m)")

figure(2)
plot(t,y)
title('Y ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Y (m)")

% Localizar m�ximos
imax = 0;
for k= 2:N-1 
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

%Mostrar resultados
fprintf(' Per�odo m�todo num�rico:   %d \n', T);
freqAng= 1/T
