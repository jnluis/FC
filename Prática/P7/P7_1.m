%% Problema 7.1 - Poço de potencial infinito a uma dimensão
clc;
clear all;
close all;

a=1; % metade da largura do poço
h=0.001;
x=-a:h:a;
N=numel(x);

n=1;
%n=2;

En=n^2*pi^2/(8*a^2);

% -------------- Numerov progressivo e normalização ------------------
% Primeiros dois valores de psi
psi1=0;
psi2=h; % Tanto faz

    g=2*En; % Neste problema, g não varia com x
    % Constantes auxiliarres para o método de Numerov

    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);

    psi=zeros(1,N);
    psi(1)=psi1;
    psi(2)=psi2;

    for n=2:N-1
        psi(n+1)=(-aux1*psi(n-1)+aux2*psi(n))/aux1;
    end

figure(1)
plot(x,psi)
xlabel('x')
ylabel('\Psi')

% Normalização
c_norm= sqrt(trapz(x,psi.^2)); % Integração numérica usando a regra dos trapézios
psi_norm=psi/c_norm;
figure(2)
plot(x,psi_norm)
xlabel('x')
ylabel('\Psi_{norm}')

figure(3)
plot(x,abs(psi_norm.^2))
xlabel('x')
ylabel('\Psi_norm^2')