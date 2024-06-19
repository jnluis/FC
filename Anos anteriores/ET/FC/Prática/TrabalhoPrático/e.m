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
%*****************************Alínea e)***********************************
clc,clear,close all

%%  Variáveis

M=1;    %massa
w=3;
lameda=0.005;
n=1;
a1=2;
b=4*a1;

%vetor x 
h=0.01;    
x=-b:h:0;       
N1=length(x);
N=2*N1;
x1=-b:h:b;

%vetores 
y=zeros(1,N);
y_pro=zeros(1,N1);

g=zeros(1,N);
%potencial
V=0.5*w^2*x.^2+lameda*x.^4;

%Condições fronteira
y(1) = 0;
y(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y(N-1) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov
y(N) = 0;

y_pro(1) = 0;
y_pro(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov

%Parametros para o shooting
B=0;        
tol=10^-10;
comparador=1;

%vetores para shooting
Ns=100;
result=zeros(1,Ns);

%Estimar os valores para o shooting
Enharm=(n+1/2)*w;   %Componente harmonica
%condiçoes para o shooting  (10.5 valor base para n=3)
E(1)=4.49;
E(2)=4.51;

%%  Método de shooting e Numerov

for k=1:Ns

    g=-2.*(V-E(k));

    C1=1+h^2/12*g; 
    C2=2*(1-5*h^2/12*g);

    %Metodo de numerov progressivo
    for i=2:N1-1
        y_pro(i+1)=C1(i+1)^(-1)*(-C1(i-1)*y_pro(i-1) + C2(i)*y_pro(i));
    end

    result(k)=y_pro(end);

    %Aplicação do metodo da secante  (E-> guess ; y_pro-> result)
    if(k>1)
        m=(result(k)-result(k-1))/(E(k)-E(k-1));

        if (m == 0)
          break;
        end

        E(k+1)=E(k)+(B-result(k))/m;      %guess(3) = guess(2) + (B-result(2))/m 

        %condição de paragem
        if(abs(E(k+1)-E(k))<tol)
          break;
        end
    end  
end       

%Normalização
c_a=trapz(x,y_pro.^2);      
ynor_pro=y_pro./sqrt(2*c_a);

%criação de ynormalizado completo
ynor_reg=-flip(ynor_pro);
ynor=[ynor_pro(1:N1),ynor_reg(2:N1)];

%%  Gráfico

fprintf('E1 = %f \n',E(end)) %Valor próprio da energia 

plot(x1,ynor)    
grid on
xlabel('x')
ylabel('\psi_{norm}')
title("\psi_{norm} em função de x")