%4)
clc;
clear all;
close all;

%variáveis
L=0.25;
C=10^-3;
a=1/(L*C);

Hhs = [10^-1 10^-2 10^-3 10^-4 ];
Nh = length(Hhs); 

for j=1:Nh
    h = Hhs(j);

        t=0:h:0.5;
    
        %equações
        N=length(t);
        q=zeros(N,1);
        aux_v=zeros(N,1);
        v=zeros(N,1);
        v(1)=5;
    
            for k=1:N-1
                aux_v(k+1)=aux_v(k) - a*v(k)*h;
                v(k+1)=v(k)+aux_v(k)*h;
            end
    
        i = aux_v * C;
        q=v*C;
    
    figure()
    subplot(3,1,1)
    plot(t,q);
    grid on;
    title('carga');
    ylabel('c/coulomb');
    xlabel('t/s');

    subplot(3,1,2)
    plot(t,i);
    grid on;
    title('corrente');
    ylabel('a/ampere');
    xlabel('t/s');


    subplot(3,1,3)
    plot(t,v);
    grid on;
    title('tensão');
    ylabel('v/volt');
    xlabel('t/s');

    vmax =  max(v);
    erro_global(j) = abs(v(1)-vmax);
end

figure()
plot(log(Hhs),log(erro_global));

figure()
scatter(log(Hhs),log(erro_global));
lsline

