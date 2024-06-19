%% Problema 2.2 - Órbita de Mercúrio
clc;
close all;
clear all;

t0= 0;
tf= 0.5;
h= 0.0001;
t = t0:h:tf;
N = numel(t);

x0 = 1.0167; % AU -> unidade astronomica
y0= 0;
vx0 = 0;
vy0 = 8.2; % AU/ano

% Constantes físicas (vêm do proprio problema)
GMs = 4*pi^2;

x = zeros(1,N);
x(1) = x0;
y = zeros(1,N);
y(1) = y0;
r = zeros(1,N);
r(1) = norm([x0 y0]);
vx = zeros(1,N);
vx(1) = vx0;
vy = zeros(1,N);
vy(1) = vy0;
v = zeros(1,N);
v(1) = norm([vx0 vy0]);
ang =  zeros(1,N);
ang(1) = mod(atan2(y0,x0), 2*pi); % = 0

ax = zeros(1,N);
ax(1) = 0;
ay = zeros(1,N);
ay(1) = 0;

for k=1:N-1   
    vx(k+1) = vx(k) -GMs*x(k)/r(k)^3 *h;
    vy(k+1) = vy(k) -GMs*y(k)/r(k)^3 *h;
    v(k+1) = norm([vx(k+1) vy(k+1)]); % em x e y


    x(k+1) = x(k) + vx(k+1) *h;
    y(k+1) = y(k) + vy(k+1) *h;
    r(k+1) = norm([x(k+1) y(k+1)]); % em x e y

    ang(k+1) = mod(atan2(y(k+1),x(k+1)), 2*pi);
    
    if ang(k+1) < ang(k)
        break
    end
end

% Eliminar os zeros finais
N = k+1;
t = t(1:N);
x = x(1:N);
y = y(1:N);
r = r(1:N);
vx = vx(1:N);
vy = vy(1:N);
v = v(1:N);
ang = ang(1:N);
ang(N) = ang(N) + 2* pi;  % Para as contas darem certo na alinea b),
                            % o angulo não pode dar um salto de -2 pi

figure(1)
plot(x,y)
hold on
scatter(0,0,30, 'r', 'filled')
title('Alinea a ')
%axis([-1.5 1.5 -1.5 1.5])
xlabel('x')
ylabel('y')
set(gca,'PlotBoxAspectRatio',[1 1 1])

% b)
T = interp1(ang(end-1:end), t(end-1:end), 2*pi);
fprintf('Periodo = %d\n\n', T)