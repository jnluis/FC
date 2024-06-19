%exercicio 3)
clc;
clear all;
close all;

%variáveis
tfin=0.5;
L=0.25;
C=10e-3;
a = 1/(L*C);
vo=5;

%valor teórico
w = 1/sqrt(L*C);
vc_tfin_a = vo*cos(w.*tfin);

%valores de h
stability = 2*sqrt(2)/w;
%os valores de h tem de ser menores que o valor de stability
h_vec=[0.001 0.002 0.005 0.02 0.01];
H=length(h_vec);

for j=1:H
    h=h_vec(j);
    t = 0:h:tfin;
    N =  length(t);
    aux = zeros(1,N);
    vc = zeros(1,N);
    vc(1) = vo;
    for  k=1:N-1
        aux(k+1) = aux(k) -a*vc(k)*h;
        vc(k+1) = vc(k) + aux(k+1)*h;
    end
    erro(j)= abs(vc(N) - vc_tfin_a);
end

aux=polyfit(log(h_vec),log(erro),1);
declive=aux(1)

figure()
plot(log(h_vec),log(erro),'.r','MarkerSize',20);
xlabel('ln(h)');
ylabel('ln(erro)');
lsl=lsline;
lsl.Color = 'k';
lsl.LineWidth = 1.2;

