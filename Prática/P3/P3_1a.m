%% Problema 3.1 - Sistema Massa/Mola
clc
clear all 
close all


x0 = 1;
v0 = 0;
k = 16;
m = 1;
t0 = 0;
tf = 15;
h = 0.01;
t = t0:h:tf;
N = length(t);
v = zeros(1, N);
x = zeros(1, N);
x(1) = x0;
v(1) = v0;


% fx = @(t,X,V) V; 
% fv = @(t,X,V) -k*X/m; 
fx = @(V) V;
fv = @(X) -k*X/m;


for i=1:N-1
    
%     r1v = -(k/m)*x(i);
%     r1x = v(i);
%     r2v = -(k/m)*(x(i)+r1x*(h/2));
%     r2x = v(i)+r1v*(h/2);
% 
%     r1x = fx(t(i),x(i),v(i));
%     r1v = fv(t(i),x(i),v(i));
%     r2x = fx(t(i)+h/2,x(i)+r1x*h/2,v(i)+r1v*h/2);
%     r2v = fv(t(i)+h/2,x(i)+r1x*h/2,v(i)+r1v*h/2);
    
    r1x = fx(v(i));
    r1v = fv(x(i));
    r2x = fx(v(i)+r1v*h/2);
    r2v = fv(x(i)+r1x*h/2);
    
    v(i+1) = v(i)+r2v*h;
    x(i+1) = x(i)+r2x*h;
    
end

xa = x0*cos(sqrt(k/m)*t);
va = -x0*sqrt(k/m)*sin(sqrt(k/m)*t);

Ec = (1/2)*m*v.^2;
Ep = (1/2)*k*x.^2;
Et = Ec + Ep;
Eca = (1/2)*m*va.^2;
Epa = (1/2)*k*xa.^2;
Eta = Eca + Epa;

figure(1)
plot(t, x, 'b-')
hold on 
plot(t, xa, 'k-')
hold off

figure(2)
plot(t, v, 'r-')
hold on 
plot(t, va, 'k-')
hold off 

figure(3)
plot(x, v, 'm-')
hold on
plot(xa, va, 'k-')
hold off

figure(4)
plot(t, Et, 'g-')
hold on 
plot(t, Eta, 'k-')
hold off




