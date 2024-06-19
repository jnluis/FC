%problema 3.1
clc
close all
clear all

v_inicial = 0;%m/s
K = 16;%N/m
x_inicial = 1;%m
m = 1;%Kg

%analitico

w = sqrt(K/m);
h = 0.01;

t_analitico = 0:h:10;
Namostras = length(t_analitico);
x_analitico= x_inicial * cos(w*t_analitico);

figure()
plot(t_analitico,x_analitico)
legend('analitico')
ylabel('distância(m)')
xlabel('tempo(s)')

%
t = 0:h:10;
Namostras = length(t);
x = zeros(Namostras,1);
x(1)=x_inicial;

v = zeros(Namostras,1);
v(1)=v_inicial;

for k=1:Namostras-1

    r1_v = -K/m*x(k); 
    r1_x = v(k);
    r2_v = -K/m*(x(k) + 0.5*r1_x*h) ;
    r2_x = v(k) + 0.5*r1_v*h;

    v(k+1) = v(k) + r2_v*h;
    x(k+1) = x(k) + r2_x*h;

end

x(end)

figure()
plot(t,x)
legend('metodo')
ylabel('distância(m)')
xlabel('tempo(s)')
%%
%alínea c
clc
close all
clear all

v_inicial = 0;%m/s
K = 16;%N/m
x_inicial = 1;%m
m = 1;%Kg

fv = @(T,X,V)  -K*X/m ;
fx = @(T,X,V)   V;

w = sqrt(K/m);
h = 0.01;

t = 0:h:10;
Namostras = length(t);
x = zeros(Namostras,1);
x(1)=x_inicial;

v = zeros(Namostras,1);
v(1)=v_inicial;

for k=1:Namostras-1

    r1_v = fv(t(k),x(k),v(k));
    r1_x = fx(t(k),x(k),v(k));
    r2_v = fv(t(k)+h/2,x(k)+ 0.5*r1_x*h,v(k)+ 0.5*r1_x*h);
    r2_x = fx(t(k)+h/2,x(k)+ 0.5*r1_x*h,v(k)+ 0.5*r1_v*h);

    v(k+1) = v(k) + r2_v*h;
    x(k+1) = x(k) + r2_x*h;

end

figure()
plot(t,x)
legend('metodo')
ylabel('distância(m)')
xlabel('tempo(s)')
