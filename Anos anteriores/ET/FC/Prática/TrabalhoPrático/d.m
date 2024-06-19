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
%*****************************Alínea d)***********************************
clc,clear,close all

%%  Variáveis

lameda=0.005;
w=3;

n=[0,1,2,3,4];

%valores de dE obtidos da alínea c)
dE=[4.1639e-4 , 0.0021 , 0.0054 , 0.0104 , 0.0170];

%ajuste polinomial (a*x^2+bx+c)
p=polyfit(n,dE,2);
a=p(1);
b=p(2);
c=p(3);

%pontos para gráfico do ajuste polinomial
x=linspace(0,4,100);
dEp=a*x.^2+b*x+c;

%dE pela Teoria das perturbações
dEteo=3/4*lameda/w^2*(2*n.^2+2*n+1);

%%  Gráfico

figure(1)
plot(n,dE,'ro',x, dEp, 'b-'); 
xlabel('n');
ylabel('dE');
title('Ajuste Polinomial de Segunda Ordem e valores de dE obtidos em c)');
legend('Valores obtidos em c)', 'Ajuste Polinomial');
grid on;

figure(2)
plot(n,dE,'ro',n,dEteo,'b*',x,dEp,'m');
xlabel('n');
ylabel('dE');
title('Comparação entre os valores obtidos em c) e os obtidos pela teoria das perturbações');
legend('Valores Teoria das perturbações', 'Valores obtidos em c)','Ajuste Polinomial');
grid on;
