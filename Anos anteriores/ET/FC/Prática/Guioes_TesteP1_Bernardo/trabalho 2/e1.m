clc;close all; clear;
%% A. Oscilador Harmónico Simples - Euler e Euler-Cromer

K  =  1; %N/m
m  = 1;  %Kg
x0 = 1;  %m

%c)
t0 = 0; tf = 50; h = 0.01; %s
t = t0:h:tf;
N = length(t);

v = zeros(1,N);
x = zeros(1,N);

v1 = zeros(1,N);
x1 = zeros(1,N);

x(1) = x0;
v(1) = 0;
x1(1) = x0;
v1(1) = 0;

for k = 1:N-1
    v(k+1)=v(k) + (-K/m*x(k)).*h;
    x(k+1)= x(k) + v(k).*h;
end

for k = 1:N-1
    v1(k+1)= v1(k) + (-K/m*x1(k)).*h;
    x1(k+1)= x1(k) + v1(k+1).*h; 
end

w = sqrt(K/m);
xx = x0*cos(w.*t);
vx = -x0*w*sin(w.*t);


subplot(2,1,1)
plot(t,v,t,v1,'--r',t,vx,':k')
title('Velocidade');xlabel('t / s');ylabel('v / (m/s)')
legend('M. Euler','M. Euler-Cromer','Sol. Analítica')

subplot(2,1,2)
plot(t,x,t,x1,'--r',t,xx,':k')
title('Posição');xlabel('t / s');ylabel('x / m')
legend('M. Euler','M. Euler-Cromer','Sol. Analítica')


%% Euler Implicito
v2 = zeros(1,N);
x2 = zeros(1,N);

x2(1) = x0;
v2(1) = 0;

for k = 1:N-1
    x2(k+1)= (x2(k) + v2(k)*h)/(1+w^2*h^2);
    v2(k+1)= v2(k)+(-K/m*x2(k+1))*h;
end

figure(3)

subplot(2,1,1)
plot(t,v,t,v2,'--r',t,vx,':k')
title('Velocidade');xlabel('t / s');ylabel('v / (m/s)')
legend('M. Euler','M. Euler Implicito','Sol. Analítica')

subplot(2,1,2)
plot(t,x,t,x2,'--r',t,xx,':k')
title('Posição');xlabel('t / s');ylabel('x / m')
legend('M. Euler','M. Euler Implicito','Sol. Analítica')

%% Crank-Nicolson
v3 = zeros(1,N);
x3 = zeros(1,N);

x3(1) = x0;
v3(1) = 0;

A = [1 h/2*w^2; -h/2 1];

for k=1:N-1
    B = [v3(k)-h/2*w^2*x3(k);1/2*v3(k)*h+x3(k)];
    LU = linsolve(A,B);
    v3(k+1) = LU(1);
    x3(k+1) = LU(2);
end
figure(4)

subplot(2,1,1)
plot(t,v,t,v3,'--r',t,vx,':k')
title('Velocidade');xlabel('t / s');ylabel('v / (m/s)')
legend('M. Euler','M. Crank-Nicolson','Sol. Analítica')

subplot(2,1,2)
plot(t,x,t,x3,'--r',t,xx,':k')
title('Posição');xlabel('t / s');ylabel('x / m')
legend('M. Euler','M. Crank-Nicolson','Sol. Analítica')


%% Energia Mecanica
figure(2)
Em_a = 1/2*K*x0;
Em_1 = 1/2*K.*x1.^2 + 1/2*m.*v1.^2;
Em_2 = 1/2*K.*x.^2 + 1/2*m.*v.^2;
Em_3 = 1/2*K.*x2.^2 + 1/2*m.*v2.^2;
Em_4 = 1/2*K.*x3.^2 + 1/2*m.*v3.^2;
subplot(4,1,2)
plot(t,Em_1)
title('Energia Mecânica (Euler-Cromer)');xlabel('t / s');ylabel('Em / J')
subplot(4,1,1)
plot(t,Em_2)
title('Energia Mecânica(Euler)');xlabel('t / s');ylabel('Em / J')
subplot(4,1,3)
plot(t,Em_3)
title('Energia Mecânica (Euler Implicito)');xlabel('t / s');ylabel('Em / J')
subplot(4,1,4)
plot(t,Em_4)
title('Energia Mecânica (Euler Implicito)');xlabel('t / s');ylabel('Em / J')

%% alinea e) com Euler-Cromer
xmax = [];
tmax = [];
imax = 0;

for k = 2:N-1
    if and(x1(k+1)-x1(k)<=0,x1(k)-x1(k-1)>=0) % condição de máximo
        imax = imax+1;
        POL = lagr(t(k-1:k+1),x1(k-1:k+1));
        tmax(imax) = POL(1);
        xmax(imax) = POL(2);
    end
end

A = mean(xmax);
NI = length(tmax);

PLF = polyfit(1:NI,tmax,1);
figure()
plot(1:NI,tmax)
T = PLF(1)



