%2)
%d)
clc;
clear all;
close all;

%variáveis
L=0.25;
C=10^-3;
a=1/(L*C);

%tempo
h=0.001;
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
%%
%e)
clc;
clear all;
close all;

%variáveis
L=0.25;
C=10^-3;
a=1/(L*C);

%tempo
h=0.001;
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

%valores teóricos
v0 = 5;
q0 = v0*C;
w = 1/sqrt(L*C);

v_a = v0*cos(w.*t);
q_a = q0*cos(w.*t);
i_a = -q0 * w * sin(w.*t);

figure()
subplot(3,1,1)
plot(t,q,t,q_a);
grid on;
title('carga');
ylabel('c/coulomb');
xlabel('t/s');
legend('numérico','analítico');

subplot(3,1,2)
plot(t,i,t,i_a);
grid on;
title('corrente');
ylabel('a/ampere');
xlabel('t/s');
legend('numérico','analítico');


subplot(3,1,3)
plot(t,v,t,v_a);
grid on;
title('tensão');
ylabel('v/volt');
xlabel('t/s');
legend('numérico','analítico');
%%
%f)
clc;
clear all;
close all;

%variáveis
L=0.25;
C=10^-3;
a=1/(L*C);

%tempo
h=0.001;
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

%tempo
h=0.0001;
t2=0:h:0.5;

%equações
N=length(t2);
q2=zeros(N,1);
aux_v2=zeros(N,1);
v2=zeros(N,1);
v2(1)=5;

for k=1:N-1
    aux_v2(k+1)=aux_v2(k) - a*v2(k)*h;
    v2(k+1)=v2(k)+aux_v2(k)*h;
end

i2 = aux_v2 * C;
q2=v2*C;

figure()
subplot(3,2,1)
plot(t,q);
grid on;
title('carga');
ylabel('c/coulomb');
xlabel('t/s');
title('0.001');

subplot(3,2,3)
plot(t,i);
grid on;
title('corrente');
ylabel('a/ampere');
xlabel('t/s');


subplot(3,2,5)
plot(t,v);
grid on;
title('tensão');
ylabel('v/volt');
xlabel('t/s');

subplot(3,2,2)
plot(t2,q2);
grid on;
title('carga');
ylabel('c/coulomb');
xlabel('t/s');
title('0.0001');

subplot(3,2,4)
plot(t2,i2);
grid on;
title('corrente');
ylabel('a/ampere');
xlabel('t/s');


subplot(3,2,6)
plot(t2,v2);
grid on;
title('tensão');
ylabel('v/volt');
xlabel('t/s');
%%
%g)
clc;
clear all;
close all;

%variáveis
L=0.25;
C=10^-3;
a=1/(L*C);

%tempo
h=0.001;
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

%valores maximos e respetivo t
periodo = islocalmax(v);
periodo = t(periodo);
periodo = diff(periodo);
%periodo é aproximadamente 0.1s

%para um grafico bonito
indices = islocalmax(v);

figure()
plot(t,v,t(indices),v(indices),'r*');
grid on;
title('tensão');
ylabel('v/volt');
xlabel('t/s');
%%
%h)
clc;
clear all;
close all;

%variáveis
L=0.25;
C=10^-3;
a=1/(L*C);

%tempo
h=0.001;
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

q0 = q(1);
e_condensador = 0.5* ((q.^2)/C);
e_bobina = 0.5 * L * i.^2;
e_sum = 0.5 * ((q0.^2/C));
e_total = e_bobina + e_condensador;

figure()
plot(t,e_total)