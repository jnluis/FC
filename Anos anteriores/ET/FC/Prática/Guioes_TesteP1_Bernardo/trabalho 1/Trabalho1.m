%13/02/2023

 clc;
 clear all;
 close all;
% %--------------------dá merda com z a apontar para baixo--------------------
% %v - velocidade
% %z - posição
% %t - tempo
% 
% %vi = 0 (velocidade inicial)
% 
% m = 150*10^-3; %kg
% g = 9.8;
% a = 3;
% 
% %a) b) c) nos apontamentos
% 
% %d)
% 
% ti=0;
% tf=1.5;
% h=0.2;
% t=ti:h:tf;
% N=length(t);
% 
% v=zeros(1,N);
% z=zeros(1,N);
% 
% %condições iniciais
% zi = 0;
% z(1)=0;%introduzir posição inicial 
% vi = 0; %introduzir velocidade inicial 
% v(1)=vi;
% 
% for k=1:N-1
%     v(k+1) = v(k) + g*h;
%     z(k+1) = z(k) + v(k)*h;
%     if(z(k+1) >6)
%         break;
%     end
% end
% 
% figure()
% plot(t,v);
% title('Gráfico da velocidade em função do tempo(c/z a apontar para baixo)');
% xlabel('t/s');
% ylabel('v/ms^{-1}');
% 
% figure()
% plot(t,z);
% title('Gráfico da position em função do tempo(c/z a apontar para baixo)');
% xlabel('t/s');
% ylabel('z/m');
% 
% t_solo=interp1(z(k:k+1),t(k:k+1),0,'linear');
% v_solo=interp1(z(k:k+1),v(k:k+1),0,'linear');
% t(k+1)=t_solo;
% z(k+1)=6;
% v(k+1)=v_solo;

%%
%------------------ *com z a apontar para cima* --------------------------
% clc;
% clear all;
% close all;

m = 150*10^-3; %kg
g = -9.8;
a = 3;

%a) b) c) nos apontamentos

%d)

ti=0;
tf=1.5;
h=0.2;
t=ti:h:tf;
N=length(t);

v=zeros(1,N);
z=zeros(1,N);

%condições iniciais
zi = 6;
z(1)=zi;%introduzir posição inicial 
vi = 0; %introduzir velocidade inicial 
v(1)=vi;

for k=1:N-1
    v(k+1) = v(k) + g*h;
    z(k+1) = z(k) + v(k)*h;
    if(z(k+1) < 0)
        break;
    end
end

figure()
plot(t,v);
title('Gráfico da velocidade em função do tempo c/z a apontar para cima');
xlabel('t/s');
ylabel('v/ms^{-1}');

figure()
plot(t,z);
title('Gráfico da position em função do tempo(c/z a apontar para arriba)');
xlabel('t/s');
ylabel('z/m');

t_solo=interp1(z(k:k+1),t(k:k+1),0,'linear');
v_solo=interp1(z(k:k+1),v(k:k+1),0,'linear');
t(k+1)=t_solo;
z(k+1)=0;
v(k+1)=v_solo;

 %%
 % *problema 1.2*

clc;
close all
clear;

L = 0.25;
C = 1*10^-3;
a = 1/(L*C);

ti=0;
tf=0.5;
h=0.001;
t=ti:h:tf;
N=length(t);

v=zeros(1,N);
v(1)=5; %tensão inicial

aux_v = zeros(1,N);
i=zeros(1,N);
q=zeros(1,N);

for k=1:N-1

    aux_v(k+1) = aux_v(k) - a*v(k)*h;
    v(k+1) = v(k) + aux_v(k)*h;
  
    q(k) = C * v(k);
    i(k) = C * aux_v(k);

end

% figure()
% plot(t,aux_v,t,v);
% 
% figure()
% plot(t,q,t,i);
% title('Gráfico da tensão e corrente , não me perguntem qual é qual');

figure()
plot(t,v,t,q,t,i);
legend({'tensão','carga','corrente'},'Location','northwest');
grid on;

figure()
subplot(3,1,1);
plot(t,v);
title('gráfico da tensão')
grid on;

subplot(3,1,2)
plot(t,i);
title('gráfico da corrente')
grid on;

subplot(3,1,3)
plot(t,q);
title('gráfico da carga');
grid on;

%%
%alínea e)

 Vo = 5 ;
 Qo = Vo * C ; 
 w = 1/sqrt(L*C);

 for k=1:N

    v_teorico(k) = Vo * cos(w*t(k));
    q_teorico(k) = Qo * cos(w*t(k));
    i_teorico(k) = -Qo * w * sin(w*t(k));

 end

 figure()
subplot(3,1,1);
plot(t,v,t,v_teorico);
title('gráfico da tensão')
grid on;

subplot(3,1,2)
plot(t,i,t,i_teorico);
title('gráfico da corrente')
grid on;

subplot(3,1,3)
plot(t,q,t,q_teorico);
title('gráfico da carga');
grid on;

%%
%alínea f)

clc;
%close all
clear;

L = 0.25;
C = 1*10^-3;
a = 1/(L*C);

ti=0;
tf=0.5;
h=0.0001;
t=ti:h:tf;
N=length(t);

v=zeros(1,N);
v(1)=5; %tensão inicial

aux_v = zeros(1,N);
i=zeros(1,N);
q=zeros(1,N);

for k=1:N-1

    aux_v(k+1) = aux_v(k) - a*v(k)*h;
    v(k+1) = v(k) + aux_v(k)*h;
  
    q(k) = C * v(k);
    i(k) = C * aux_v(k);

end

% figure()
% plot(t,aux_v,t,v);
% 
% figure()
% plot(t,q,t,i);
% title('Gráfico da tensão e corrente , não me perguntem qual é qual');

figure()
plot(t,v,t,q,t,i);
legend({'tensão','carga','corrente'},'Location','northwest');
grid on;

figure()
subplot(3,1,1);
plot(t,v);
title('gráfico da tensão')
grid on;

subplot(3,1,2)
plot(t,i);
title('gráfico da corrente')
grid on;

subplot(3,1,3)
plot(t,q);
title('gráfico da carga');
grid on;

%alínea e) de f)
 Vo = 5 ;
 Qo = Vo * C ; 
 w = 1/sqrt(L*C);

 for k=1:N

    v_teorico(k) = Vo * cos(w*t(k));
    q_teorico(k) = Qo * cos(w*t(k));
    i_teorico(k) = -Qo * w * sin(w*t(k));

 end

figure()
subplot(3,1,1);
plot(t,v,t,v_teorico);
title('gráfico da tensão')
grid on;

subplot(3,1,2)
plot(t,i,t,i_teorico);
title('gráfico da corrente')
grid on;

subplot(3,1,3)
plot(t,q,t,q_teorico);
title('gráfico da carga');
grid on;

%%
%alínea g) 

indmax = find(islocalmax(v));
v_max = v(indmax);
tmax = t(indmax);

Nmax = numel(tmax);
T = (tmax(Nmax)-tmax(1))/(Nmax-1);

%c/interpolação
for imax = 1 : Nmax

    ii = indmax(imax);
    tmax_i(imax) = interp1(i(ii-1:ii+1),t(ii-1:ii+1),0);
    vmax_i(imax) = interp1(i(ii-1:ii+1),v(ii-1:ii+1),0);
end

Tint = (tmax_i(Nmax)-tmax_i(1))/(Nmax-1);
T_analitico = 2*pi/sqrt(a);

%v maximo -> corrente exatamente zero

format long

disp('merdas com interpolação')
display(tmax_i);
display(Tint);
display(T_analitico);
display(vmax_i);

format short

%%
%alinea h)
for k=1 : length(t)
    energia_condensador(k) = 0.5 * ((q(k)^2)/C);
    energia_bobina(k) = 0.5 * (L * (i(k)^2)); 
end

energia_condensador_sum = sum(energia_condensador);
energia_bobina_sum = sum(energia_bobina);

sum_energias = energia_bobina(1) + energia_condensador(1);

sum_teorico = 0.5 * ((Qo^2)/C);

%%
% * problema 1.3 *

%% 
% * problema 1.4 *


