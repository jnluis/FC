%%% Problema FR2.2: BVP – Método do Shooting
% Feito com Euler-Cromer em vez de ode45
clc;
clear all;
close all;

h= 0.001;
tf=50;
t = 0:h:tf;
N=length(t);

km=2;
m=1.5;
alpha=-0.2;
x0=1.9;
v0=0;

x=zeros(N,1);
v=zeros(N,1);

x(1)=x0;
v(1)=v0;

Tolerancia=10^-4; % a mesma tolerancia do Ex. anterior

B= -1.5; % amplitude negativa final
guess= alpha;
guess(1)=-0.19;
guess(2)=-0.2;

for is=1:150
    alpha=guess(is);
    v=zeros(N,1);
    x=zeros(N,1);

    x(1)=x0;
    v(1)=v0;

    for k=1:N-1
        % Euler-Cromer
        v(k+1)=v(k)+( (-km/m) * x(k)*(1+(3/2)*alpha*x(k)) ) *h;
        x(k+1)=x(k)+v(k+1)*h;

        im=0;
        clear xmin;
        for j=2:N-1
            if (x(j+1)> x(j) && x(j-1)>= x(j))
                im=im+1;
                aux=lagr(t(j-1:j+1),-x(j-1:j+1)); % Colocou-se o menos porque tamos
                % a ir buscar mínimos e o lagr encontra maximos
                xmin(im)=-aux(2);
            end
        end
    end
        
    AmplitudeN= mean(xmin);
    result(is)=AmplitudeN;

    if(is > 1)
        declive=(result(is)-result(is-1))/(guess(is)-guess(is-1));
        guess(is+1)= guess(is)+ (B-result(is))/declive;

        if (abs(B-result(is)) < Tolerancia)
            break;
        end
    end
end


figure(1)
plot(t,x)
xlabel('t');
ylabel('x');