clear all

%global constants
K = 1; %N/m
m = 1; %Kg
alfa = -0.1; %m^2;

options = optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

%vectors
tf = 10;
h1 = 0.01; %define hmin
h2 = 0.0099; %define hmax
ni = round(tf/h1+1);
nf = round(tf/h2+1);
N = ni:nf;

t_max = 0:h2:tf;
vx = zeros(1,length(t_max));
x = zeros(1,length(t_max));
    
h_vec = zeros(1,length(N));
xn_vec = zeros(1,length(N));

for n = 1:length(N)
    h=tf/(N(n)-1);
    t = 0:h:tf;
    
    %physical initial conditions
    x(1) = 1; %m
    vx(1) = 1; %m/s
    
    const = [h/2, K*h/(2*m), 2*alfa];
    
    for i = 1:(length(t)-1)
        func = @(xv) fcr(xv,x(i),vx(i),const);
        xv0 = [x(i),vx(i)];
        aux = fsolve(func,xv0,options);
        x(i+1) = aux(1);
        vx(i+1) = aux(2);
        
    end
    
    h_vec(n) = h;
    xn_vec(n) = x(i+1);
end
%% Plots

h_vec2 = h_vec.^2;
plot(h_vec2,xn_vec)
p = polyfit(h_vec2,xn_vec,1)