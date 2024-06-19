%% Ex2.1.
clearvars, clc ,close all;

% Input que pode ser modificado
x0 = 1;
vx0 = 0;

% Parâmetros dos métodos numéricos
h = input("Introduza o número de h em segundos: ");
tfin = 100;

% Constantes do problema
K = 1;
m = 1;

% Cálculos aux.
w = sqrt(K/m);
w2 = K/m;

% Iniciar Variaveis
t = 0:h:tfin;
N = numel(t);
x = zeros(N,1);
x(1) = x0;
vx = zeros(N,1);
vx(1) = vx0;

% Escolha do método

met = ["Euler","Euler-Cromer","Euler implícito","Crank-Nickolson"];
promptn = sprintf(['De entre os métodos,\n 1-%s,\n 2-%s,\n 3-%s,\n 4-%s,\n escolha um (1 a 4): '],met);
n = input(promptn);

switch n
    case 1
        for k = 1:N-1
            a = -K*x(k)/m;
            vx(k+1) = vx(k) + a*h;
            x(k+1) = x(k) +vx(k)*h;
        end

    case 2
        for k = 1:N-1
            a = -K*x(k)/m;
            vx(k+1) = vx(k) + a*h;
            x(k+1) = x(k) + vx(k+1)*h;
        end
    
    case 3
        A = [1 -h; w^2*h 1];
        for k = 1:N-1
            b = [x(k), vx(k)];
            [x(k+1), vx(k+1)] = linsolve(A,b);
        end
    
    case 4
         A = [1 -h/2; w^2*h/2 1];
        for k = 1:N-1
            b = [x(k)+vx(k)*h/2, vx(k)-w^2*x(k)*h/2];
            [x(k+1), vx(k+1)] = linsolve(A,b);
        end
    otherwise
        fprintf('Nope');
        return
end

%% Gráficos

%% Localizar máximos
imax = 0;
for k = 2:N-1
    if and(x(k+1)-x(k)<=0,x(k)-x(k-1)>=0)
        imax = imax+1;
        aux = lagr(t(k-1:k+1),x(k-1:k+1));
        tmax(imax) = aux(1);
        xmax(imax) = aux(2);
    end
end

% Desnecessário
nmax = imax;
plf = polyfit(1:nmax,tmax,1);
T = plf(1);
A = mean(xmax);

% Solução analitica

