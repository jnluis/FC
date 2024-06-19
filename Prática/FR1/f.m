function derivadas = f(t, solucao, m,K)
w= sqrt(K/m);
alpha = -0.1;

derivadas = zeros(2,1);
x=solucao(1);
v=solucao(2);

derivadas(1)= v;
derivadas(2)= -K*(x+ 2*alpha*x^3)/m;