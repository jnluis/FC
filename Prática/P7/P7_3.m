%% Problema 7.3 - Poço de potencial finito a uma dimensão
clc;
clear all;
close all;


% Mudar para encontrar outros valores prórpios
E(1)=1.6;
E(2)=1.8;

% Potencial
V0=20;

a=1; % metade da largura do poço
b=4; % Quando nos estendemos para os lados
h=0.001;
x=-b*a:h:b*a;
N=numel(x);

% Shooting and matching
tolera=1e-10;
nmaxE=100; % Número máximo de iterações
x_match=1.0;
ind_match=round(1+(x_match+b*a)/h);
x_left=x(1:ind_match);
x_right=x(ind_match:N);
N_left=numel(x_left);
N_right=numel(x_right);

% Potencial
V_left=zeros(1,N_left);
V_left(x_left<-a)=V0;
V_right=zeros(1,N_right);
V_right(x_right>a)=V0;

% Primeiros dois valores de psi
psi_extremo=0;
psi_seguinte_left=h/100000; % Tanto faz?
psi_seguinte_right=+psi_seguinte_left; % Tanto faz?

figure(1)
xlabel('x');ylabel('y')

for iE=1:nmaxE
    % Esquerda para a direita
    % Constantes auxiliares para o método de Numerov

    g=2*(E(iE)-V_left);
    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);

    psi_left=zeros(1,N_left);
    psi_left(1)=psi_extremo;
    psi_left(2)=psi_seguinte_right;

    for n=2:N_left-1
        psi_left(n+1)=(-aux1(n-1)*psi_left(n-1)+aux2(n)*psi_left(n))/aux1(n+1);
    end

% Direita para a esquerda
    clear g
    % Constantes auxiliares para o método de Numerov
    g= 2*(E(iE)-V_right); % Neste problema, k ao quadrado muda
    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);

    psi_right=zeros(1,N_right);
    psi_right(N_right)=psi_extremo;
    psi_right(N_right)=psi_seguinte_left;

     for n=N_right-1:-1:2
        psi_right(n-1)=(-aux1(n+1)*psi_right(n+1)+aux2(n)*psi_right(n))/aux1(n-1);
     end

    plot(x_left,psi_left,x_right,psi_right);
    xlabel('x')
    ylabel('\Psi')
    pause(0.5)

    % Acerta os valores
    ratio=psi_left(N_left)/psi_right(1);
    psi_right=psi_right*ratio;

    plot(x_left,psi_left,x_right,psi_right);
    xlabel('x')
    ylabel('\Psi')
    pause(0.5)

    D_left=(25/12*psi_left(N_left)-4*psi_left(N_left-1)+3*psi_left(N_left-2)...
        -4/3*psi_left(N_left-3)+1/4*psi_left(N_left-4))/h;
    D_right=(-25/12*psi_right(1)+4*psi_right(2)-3*psi_right(3)...
        +4/3*psi_right(4)-1/4*psi_right(5))/h;

    DLog_left=D_left/psi_left(N_left);
    DLog_right=D_right/psi_right(1);
    result(iE)=(DLog_left-DLog_right)/(DLog_left+DLog_right);

    if(iE>1)
        dif=(result(iE)-result(iE-1))/(E(iE)-E(iE-1));
        if dif == 0 % Sem isto pode dar erro para tolerancias pequenas
            break
        end
        E(iE+1)=E(iE)-result(iE)/dif;

        if(abs(E(iE+1)-E(iE)) < tolera)
            break
        end
    end
end

fprintf('Energia = %f Ha\n',E(end))

psi=[psi_left(1:N_left-1), psi_right];

% Normalização
c_norm= sqrt(trapz(x,psi.^2)); % Integração numérica usando a regra dos trapézios
psi_norm=psi/c_norm;
figure(2)
plot(x,psi_norm)
xlabel('{\it x}')
ylabel('\psi_{norm}')