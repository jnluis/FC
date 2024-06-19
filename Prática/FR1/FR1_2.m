%% Problema 2.3 - Oscilador Quártico - Métodos de Runge-Kutta
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
tf= 50;
h= 0.01;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

v = zeros(N,1);
v(1) = vx0;

const = [h/2, K*h/(2*m), 2*alpha];
options = optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

xE= x;
vE= v;

fx= @(V) V; % função anónima
fv = @(X) -K*(X+ 2*alpha*X^3)/m;

for k=1:N-1
    % Runge-Kutta           
    r1v = fv(x(k));
    r1x = fx(v(k));

    r2v = fv(x(k) +r1x*h/2);
    r2x = fx(v(k) +r1v*h/2);

    r3v = fv(x(k) +3*r2x*h/4);
    r3x = fx(v(k) +3*r2v*h/4);

    v(k+1) = v(k) + (2*r1v+3*r2v+4*r3v)*h/9;
    x(k+1) = x(k) + (2*r1x+3*r2x+4*r3x)*h/9;

    % Euler-Cromer
    a = -K*( xE(k) + 2* alpha* xE(k).^3 )/m;
    vE(k+1) = vE(k) +a *h;
    xE(k+1) = xE(k) + vE(k+1) *h;
end

Em = 0.5*(K*(x.^2+alpha*x.^4)+ m*v.^2);
figure(1)
subplot(1,2,1)
plot(t,x)
title('RK3- Oscilador Quártico');
xlabel("Tempo (s)")
ylabel("Posição (m)")

subplot(1,2,2)
plot(t,Em)
title('RK3- Oscilador Quártico');
xlabel("Tempo (s)")
ylabel("Em (J)")

% Energia Mecânica de Euler-Cromer
Ec = m*vE.^2;
Ep = K*xE.^2 .* (1+alpha*(xE.^2));
Em = 0.5* (Ec + Ep);

figure(2)
subplot(1,2,1)
plot(t,xE)
title('Euler-Cromer');
xlabel("Tempo (s)")
ylabel("Posição (m)")

subplot(1,2,2)
plot(t,Em)
title('Euler-Cromer');
xlabel("Tempo (s)")
ylabel("Em (J)")