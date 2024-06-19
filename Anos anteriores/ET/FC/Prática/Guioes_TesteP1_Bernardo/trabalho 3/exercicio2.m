%2)
close all;clc;clear all
%constantes
K=16;
m=1;

%valores iniciais
vo=0;
xo=1;

%variável independente
h=0.1;
t=0:h:10;
N=length(t);

%vari+aveis dependentes
x=zeros(1,N);
v=zeros(1,N);
x(1)=xo;
v(1)=vo;

%solução analítica
w  = sqrt(K/m);
xx = xo*cos(w.*t);
vx = -xo*w*sin(w.*t);

%solução numérica
%equações das derviadas
fx = @(V) V;           %derivada de x em função do tempo
fv = @(X) -K*X/m;      %derivada da v em função do tempo

for k=1:N-1
    %neste casoo podiamos só por um dos argumentos pq a formula da derivada
    %só usa um também, teria de mudar a função anónima 

    r1x = fx(v(k));
    r1v = fv(x(k));

    r2x = fx(v(k)+r1v*h/2);
    r2v = fv(x(k)+r1x*h/2);

    r3x = fx(v(k)+r2v*h/2);
    r3v = fv(x(k)+r2x*h/2);

    r4x = fx(v(k)+r3v*h);
    r4v = fv(x(k)+r3x*h);

    x(k+1)=x(k) + h*(r1x/6+ r2x/3 + r3x/3 + r4x/6);
    v(k+1)=v(k) + h*(r1v/6+ r2v/3 + r3v/3 + r4v/6);
end

%plots
figure()
subplot(1,2,1)
plot(t,xx,t,x);
xlabel('tempo(s)');
ylabel('distâcia(m)');
legend('sol.analítica','sol.numérica');
title('posição');
%escolhi um h tão bom que os valores são literalmente são iguais sou bué
%slay
subplot(1,2,2)
plot(t,vx,t,v);
xlabel('tempo(s)');
ylabel('velocidade(m/s)');
legend('sol.analítica','sol.numérica');
title('velocidade');
%%
figure()
plot(x,v);

Em_a = (0.5*K.*xo.^2 + 1/2*m.*vo.^2)*ones(N,1);
Em_rk = 1/2*K.*x.^2 + 1/2*m.*v.^2;

figure()
plot(t,Em_a,t,Em_rk);
xlabel('tempo(s)');
ylabel('energia(j)');
legend('sol.analítica','sol.numérica');
%%
%a alínea c) era para fzr isto mas com valores de h diferentes e n era
%perciso ver o grafico e etc era so msm comparar

%alinea d)
stability = 2*sqrt(2)/w;

close all;clc;clear all
%constantes
K=16;
m=1;
tfin = 10;
%valores iniciais
vo=0;
xo=1;

%variável independente
hvec=[0.01 0.02 0.05 0.2 0.1];
H=length(hvec);
  %solução analítica
    w  = sqrt(K/m);
    x_tfin_sa = xo*cos(w.*tfin);
    %vx = -xo*w*sin(w.*t);
%solução numérica
%equações das derviadas
fx = @(V) V;           %derivada de x em função do tempo
fv = @(X) -K*X/m;      %derivada da v em função do tempo

for j=1:H
    h=hvec(j);
    t = 0:h:tfin;
    N = length(t);

    %variaveis dependentes
    x=zeros(1,N);
    v=zeros(1,N);
    x(1)=xo;
    v(1)=vo;

  
    for k=1:N-1
        %neste casoo podiamos só por um dos argumentos pq a formula da derivada
        %só usa um também, teria de mudar a função anónima

        r1x = fx(v(k));
        r1v = fv(x(k));

        r2x = fx(v(k)+r1v*h/2);
        r2v = fv(x(k)+r1x*h/2);

        r3x = fx(v(k)+r2v*h/2);
        r3v = fv(x(k)+r2x*h/2);

        r4x = fx(v(k)+r3v*h);
        r4v = fv(x(k)+r3x*h);

        x(k+1)=x(k) + h.*(r1x/6+ r2x/3 + r3x/3 + r4x/6);
        v(k+1)=v(k) + h.*(r1v/6+ r2v/3 + r3v/3 + r4v/6);
    end
     erro(j)= abs(x(N) - x_tfin_sa);
%     %plots
%     figure()
%     subplot(1,2,1)
%     plot(t,xx,t,x);
%     xlabel('tempo(s)');
%     ylabel('distâcia(m)');
%     legend('sol.analítica','sol.numérica');
%     %escolhi um h tão bom que os valores são literalmente são iguais sou bué
%     %slay
%     subplot(1,2,2)
%     plot(t,vx,t,v);
%     xlabel('tempo(s)');
%     ylabel('velocidade(m/s)');
%     legend('sol.analítica','sol.numérica');
   
end

aux=polyfit(log(hvec),log(erro),1);
declive=aux(1)

figure()
plot(log(hvec),log(erro),'.r','MarkerSize',20);
xlabel('ln(h)');
ylabel('ln(erro)');
lsl=lsline;
lsl.Color = 'k';
lsl.LineWidth = 1.2;

% xL=xlim;
% yL=ylim;
% posX = xL(1) + 0.15*(xL(2)-xL(1));
% posY = yL(1) - 0.25*(yL(2)-yL(1));
% text(posx,posy,sprintf('declive= %.3f'))


