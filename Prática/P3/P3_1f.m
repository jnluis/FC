%% Problema 3.1 - Sistema Massa/Mola
clc;
close all;
clear all;

x0 = 1; % posição inicial
v0 = 0;

K  = 16; %N/m
m = 1; % 1 kg

w = sqrt(K/m);
w2 = K/m;

t0= 0;
tf= 15;
h= 0.01;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

v = zeros(N,1);
v(1) = v0;

xE= x;
vE= v;

fx= @(V) V; % função anónima
fv = @(X) -K*X/m;

for k=1:N-1
    % Runge-Kutta           
    r1v = fv(x(k));
    r1x = fx(v(k));

    r2v = fv(x(k) +r1x*h/2);
    r2x = fx(v(k) +r1v*h/2);

    v(k+1) = v(k) + r2v*h;
    x(k+1) = x(k) + r2x*h;

    % Euler
    xE(k+1) = xE(k) + vE(k)*h;
    vE(k+1) = vE(k) -K*xE(k)/m*h;
end

Ec = 0.5 * m * vE.^2;
Ep = 0.5 * K * xE.^2;
EE = Ec + Ep; % Euler

E = 0.5*K*x.^2+0.5*m*v.^2; % Runge-Kutta

% Solução analítica
Esa = 0.5*(K*x0^2+m*v0^2)*ones(N,1);
xsa= x0*cos(w*t);
vsa = w*x0*sin(w*t);

figure(1)
subplot(2,3,1)
plot_1= plot(t,x,'k',t,xsa,'r');
title('x(t)');

subplot(2,3,2)
axis off
hl = legend(plot_1,{"Runge-Kutta", "Sol. analitica"});
newPosition = [0.43 0.72 0.15 0.15];
newUnits= 'normalized';
set(hl, 'Position', newPosition,'Units', newUnits);

subplot(2,3,3)
plot(t,v,'k',t,vsa,'r')
title('v(t)')

subplot(2,2,3)
plot(v,x,'k',vsa,xsa,'r')
title('v(x)')

subplot(2,2,4)
plot(t,E,'k',t,Esa,'r')
title('E(t)')

%% Figura 2

figure(2)
subplot(2,3,1)
plot_1= plot(t,x,'k',t,xE,'r');
title('x(t)');

subplot(2,3,2)
axis off
hl = legend(plot_1,{"Euler", "Sol. analitica"});
newPosition = [0.43 0.72 0.15 0.15];
newUnits= 'normalized';
set(hl, 'Position', newPosition,'Units', newUnits);

subplot(2,3,3)
plot(t,v,'k',t,vE,'r')
title('v(t)')

subplot(2,2,3)
plot(v,x,'k',vE,xE,'r')
title('v(x)')

subplot(2,2,4)
plot(t,EE,'k',t,Esa,'r')
title('E(t)')

