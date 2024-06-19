%*************************************************************************
%
% NOME 1: Bruno Figueiredo
% MEC 1: 103489
% Turma: PL8
% 
%*************************************************************************
%
% NOME 2: Laura Villalba
% MEC 2: 102847
% Turma: PL8
% 
%*************************************************************************
%
% NOME 2: Rafael Morgado
% MEC 2: 104277
% Turma: PL8
%
%*****************************Alínea c)***********************************
clc,clear,close all

%%  Variáveis

M=1;    %massa
w=3;
lameda=0.005;
n=[0,1,2,3,4];
a1=2;
b=4*a1;

%vetor x 
h=0.01;    
x=-b:h:b;       
N=length(x);

%vetores 
y=zeros(1,N);
y_pro=zeros(1,N);
y_reg=zeros(1,N);

g=zeros(1,N);
%potencial
V=0.5*w^2*x.^2+lameda*x.^4;

x_match=round(2*N/5);       %indice do centro aprox

%Condições fronteira
y(1) = 0;
y(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y(N-1) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y(N) = 0;

y_pro(1) = 0;
y_pro(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y_pro(N-1) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y_pro(N) = 0;

y_reg(1) = 0;
y_reg(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y_reg(N-1) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y_reg(N) = 0;

%Parametros para o shooting
B=0;        
tol=10^-10;
k=1;
comparador=1;

%vetores para shooting
Ns=100;
result=zeros(1,Ns);
%E=zeros(1,Ns);

%Estimar os valores para o shooting
Enharm=(n+1/2)*w;   %Componente harmonica

%vetor para os vários ynormalizados e dE
ynor=zeros(5,N);
dE=zeros(1,5);
EE=zeros(1,5);

%%  Método de shooting e Numerov

for j=1:5

    E=[Enharm(j)-0.01 Enharm(j)+0.01];  %condiçoes para o shooting

    for k=1:Ns

        g=-2*(V-E(k));

        C1=1+h^2/12*g; 
        C2=2*(1-5*h^2/12*g);

        %Metodo de numerov progressivo
        for i=2:N-1
            y_pro(i+1)=C1(i+1)^(-1)*(-C1(i-1)*y_pro(i-1) + C2(i)*y_pro(i));
        end

        %derivadas do psi progressivo
        dy_pro = 1/h*(-25/12*y_pro(x_match)+4*y_pro(x_match+1)-3*y_pro(x_match+2)+4/3*y_pro(x_match+3)-1/4*y_pro(x_match+4));

        %Metodo de numerov regressivo
        for i=N-1:-1:2
            y_reg(i-1) = C1(i-1)^(-1)*(C2(i)*y_reg(i)-C1(i+1)*y_reg(i+1));
        end

        %ver ponto de descontinuidade
        ypro_xmatch=y_pro(x_match); 
        yreg_xmatch=y_reg(x_match);
    
        %Forçar a igualdade
        y_reg = y_reg*ypro_xmatch/yreg_xmatch;

        %derivadas do psi regressivo
        dy_reg = 1/h*(25/12*y_reg(x_match)-4*y_reg(x_match-1)+3*y_reg(x_match-2)-4/3*y_reg(x_match-3)+1/4*y_reg(x_match-4));

        %criar psi
        y(1:x_match) = y_pro(1:x_match);
        y(x_match:N) = y_reg(x_match:N); 

        %calculo do result para o shooting
        comparador= (dy_pro/ypro_xmatch - dy_reg/yreg_xmatch)/(dy_pro/ypro_xmatch + dy_reg/yreg_xmatch);
        result(k)=comparador;

        %Aplicação do metodo da secante  (E-> guess ; yin-> result)
        if(k>1)
            m=(result(k)-result(k-1))/(E(k)-E(k-1));

            if (m == 0)
                break;
            end

            E(k+1)=E(k)+(B-result(k))/m;      %guess(3) = guess(2) + (B-result(2))/m

            if (abs(E(k+1) - E(k)) < tol)
                break;
            end    
        end  
    end

    %guardar valores das energias
    dE(j)=E(end)-Enharm(j);
    EE(j)=E(end);

    %Normalização
    c_c=trapz(x,y.^2);      
    ynor(j,:)=y./sqrt(c_c);
end

%%  Gráficos

figure(1)
plot(x,ynor)    
grid on
xlabel('x')
ylabel('\psi_{norm}')
title("\psi_{norm} em função de x")
legend('E_0','E_1','E_2','E_3','E_4');

figure(2)
plot(x,V,'b','LineWidth',2)
grid on
hold on
ylim([0,15])
% Plot das linhas de energia
colors = {'r', 'b', 'm', 'k', 'g'}; % Cores para cada linha
for i = 1:5
    energia = EE(i);
    color_index = mod(i-1, numel(colors)) + 1; % Índice da cor
    line(xlim, [energia, energia], 'Color', colors{color_index}, 'LineStyle', '--');
end
xlabel('x');
ylabel('V(x)');
title('Poço de Potencial V(x) juntamente com as 5 Energias próprias');
legend('V(x)','E_n=0','E_n=1','E_n=2','E_n=3','E_n=4')