clear all

%global physical constants
M = 1.5; %Kg
K = 2; %N/m

%para teste da ODE45 -- depois não é usado, uma vez que é este parâmetro
%que se quer descobrir
%alfa = -0.2; % real deve dar -0.132890

%ODE45 constants
reltol = 3*10^(-14);
abstol_1 = 1*10^(-13);
abstol_2 = 1*10^(-13);

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

%system constants
t0 = 0; %s
x0 = 1.9; %m
v0 = 0; %m/s
tf = 10; %s

%shooting goal
B = -1.5; %m (amplitude)

%vectors
N = 10; %melhor dar logo memória a mais para o ciclo não  ficar lento
tol = 10^(-12);
alfa = zeros(1,N); %o que se quer descobrir
Aneg = zeros(1,N); %o que se está a condicionar

% guess(1) e guess(2)
alfa(1) = -0.1;
alfa(2) = -0.2;

%para entrar no while
comparator = 0;

%Método de Shooting
k = 1;
while abs(B - comparator) > tol
    [t,sol] = ode45(@(t,sol)f(t,sol,K,M,alfa(k)),[t0 tf],[x0 v0],options);
    
    x = sol(:,1);
    
    comparator = min(x);
    Aneg(k) = min(x);
    
    if k > 1 %para poder calcular result(1) e result(2) pelo ODE45 primeiro
        declive = (Aneg(k)-Aneg(k-1))/(alfa(k)-alfa(k-1));
        alfa(k+1) = alfa(k) + (B-Aneg(k))/declive;
    end
    
    k = k + 1;
end