%% Problema 4.2 - Método das diferenças finitas — determinação da frequência de vários modos normais de vibração
clc;
clear all;
close all;

T= 1000;
mu= 10^-3; % kg/m
L= 1; %m

w1 = (pi/L)*sqrt(T/mu);
% wn= n w1 = npi vezes 10^3

h=0.001;
x= 0:h:L;
N=numel(x);
n= N-2; % tamanho das matrizes e vetores

A=eye(n,n);
A=-2.*A;

for i=1:n-1
    A(i,i+1)=1; % diagonal acima da principal
    A(i+1,i)=1; % diagonal acima da principal
end

[vec,val] = eigs(A,3,'sm'); % 3 valores próprios e seus correspondentes vetores
% sm solicita os menores valores próprios
val_prop=diag(val);
omega= sqrt(-val_prop*T/mu)/h;

fprintf('omega do 1º modo = %f s^-1\n', omega(1));
fprintf('omega do 2º modo = %f s^-1\n', omega(2));
fprintf('omega do 3º modo = %f s^-1\n', omega(3));
fprintf('omega do 1º modo / w1 = %f\n', omega(1)/w1);
fprintf('omega do 2º modo / w1 = %f\n', omega(2)/w1);
fprintf('omega do 3º modo / w1 = %f\n', omega(3)/w1);

figure(1)
subplot(3,1,1)
plot(x,[0; vec(:,1);0])
legend('1º modo','Location','Best')
xlabel('x')
ylabel('y')
subplot(3,1,2)
plot(x,[0; vec(:,2);0])
legend('2º modo','Location','Best')
xlabel('x')
ylabel('y')
subplot(3,1,3)
plot(x,[0; vec(:,3);0])
legend('3º modo','Location','Best')
xlabel('x')
ylabel('y')
