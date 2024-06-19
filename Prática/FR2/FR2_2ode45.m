%%% Problema FR2.2: BVP – Método do Shooting
% Feito com ode45
clc
clear all
close all
fprintf('\n')

% Atenção: Há um pause no meio do programa!!!

%----------------------------------------------
% Input que pode ser modificado
x0=1.9;
v0=0.0;

% Parametros dos metodos numericos
h=0.001; % Isto é só para o output do ode45
tfin=50;
t=0:h:tfin; % Isto é só para o output do ode45
N=numel(t); % Isto é só para o output do ode45

% Constantes do problema
K=2; % K maiusculo!
m=1.5;
% Calculos auxiliares
w=sqrt(K/m);
w2=K/m;

% Dois alfas iniciais para o shooting
alfa(1)=-0.2;
alfa(2)=-0.19;
% Outros parametros de shooting
tolera=1e-12;
is_max=50;
Aneg_alvo=-1.5;

for is=1:is_max
    
    options=odeset('RelTol',1e-10,'AbsTol',[1e-10 1e-10]);
    [t,solucao]=ode45(@f_FR2_2,t,[x0 v0],options,K,m,alfa(is));
    
    x = solucao(:,1);
    v = solucao(:,2);
    
    plot(t,x)
    xlabel('{\it t}');ylabel('{\it x}')
    pause(0.5) % Para ir vendo como nos aproximamos da solucao.
    
    % Localiza minimos e calcula Aneg
    imin=0;
    clear xmin % E' preciso para nao ficarem mins de um is para outro
    for k=2:N-1
        if and(x(k+1)>x(k),x(k-1)>=x(k))
            imin=imin+1;
            aux=lagr(t(k-1:k+1),-x(k-1:k+1));
            xmin(imin)=-aux(2);
        end
    end
    Aneg(is)=mean(xmin);
        fprintf('Aneg da iteracao = %f\n',  Aneg(is))
    
    % Metodo da secante
    if(is>1)
        if abs(Aneg_alvo-Aneg(is)) < tolera
            break
        end
        dif=(Aneg(is)-Aneg(is-1))/(alfa(is)-alfa(is-1));
        alfa(is+1)=alfa(is)+(Aneg_alvo-Aneg(is))/dif;
    end
    
end

fprintf('\nalfa final = %f\n\n',  alfa(end))
