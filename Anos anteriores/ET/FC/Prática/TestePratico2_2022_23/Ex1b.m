%%% 2º Teste -> 22-23
% Feito com Euler-Cromer
clc;
clear all;
close all;

h= 0.01;
tf=100;
t = 0:h:tf;
N=length(t);

r=0.2;
omega=1;
y0=0;
Dy0=0.2;

y=zeros(N,1);
Dy=zeros(N,1);

y(1)=y0;
Dy(1)=Dy0;

Tolerancia=10^-4;
alpha=0.4;

B= -0.2; % media do valor minimo da velocidade
guess= alpha;
guess(1)=0.39;
guess(2)=0.4;

for is=1:150
    alpha=guess(is);
    Dy=zeros(N,1);
    y=zeros(N,1);

    y(1)=y0;
    Dy(1)=Dy0;

    for k=1:N-1
        % Euler-Cromer
        Dy(k+1)=Dy(k)+ ( -(y(k).^3-y(k)) + (-alpha*Dy(k)) + (r*cos(omega*t(k))) )*h;
        y(k+1)=y(k)+Dy(k+1)*h;

        im=0;
        for j=2:N-1
            if (Dy(j+1)> Dy(j) && Dy(j-1)>= Dy(j))
                im=im+1;
                aux=lagr(t(j-1:j+1),-Dy(j-1:j+1)); % Colocou-se o menos porque tamos
                % a ir buscar mínimos e o lagr encontra maximos
                ymin(im)=-aux(2);
            end
        end
    end
        
    VelocN= mean(ymin);
    result(is)=VelocN;

    if(is > 1)
        declive=(result(is)-result(is-1))/(guess(is)-guess(is-1));
        guess(is+1)= guess(is)+ (B-result(is))/declive;

        if (abs(B-result(is)) < Tolerancia)
            break;
        end
    end
end


figure(1)
plot(t,y)
xlabel('t');
ylabel('y');

figure(2)
plot(t,Dy)
xlabel('t');
ylabel('Dy');

figure(3)
plot(y,Dy)
xlabel('y');
ylabel('Dy');