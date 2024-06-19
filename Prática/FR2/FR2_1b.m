%%% Problema FR2.1: BVP – Método do Shooting
% Feito com Euler-Cromer
clc;
clear all;
close all;

h= 0.01;
tf=7;
x = 0:h:tf;
N=length(x);

beta=18;


v=1.7187;
DT0=10^-4*v;

T(1)=10^-4;
DT(1)=DT0;

Tolerancia=10^-4;

guess=v;
B= 0; % CF da derivada de T em função de x ser 0 em x=7 
guess(1)=1.9;
guess(2)=2.0;

for is=1:150
    v=guess(is);
    
    T=zeros(N,1);
    DT=zeros(N,1);
    T(1)=10^-4;
    DT(1)=10^-4*v;

    for k=1:N-1
        % Euler-Cromer
        DT(k+1) = DT(k) + (v*DT(k)-beta * exp(-1/T(k))*(1+DT(k)-v*T(k)))*h;
        T(k+1) = T(k)+DT(k+1)*h;
    end
    result(is)=DT(end); % CF

    if(is > 1)
        declive=(result(is)-result(is-1))/(guess(is)-guess(is-1));
        guess(is+1)= guess(is)+ (B-result(is))/declive;

        if (abs(B-result(is)) < Tolerancia)
            break;
        end

    end

end

figure(1)
subplot(1,2,1)
plot(x,T)
xlabel('x');
ylabel('Temperatura T');

subplot(1,2,2)
plot(x,DT)
xlabel('x');
ylabel('Derivada de T');