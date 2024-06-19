function derivadas = f(t,solucao)
% A função f retorna as derivadas dos valores do tempo e da solução da
% função f
derivadas=zeros(2,1,1); % alocar as derivadas
x=solucao(1);
y=solucao(2);
z=solucao(3);

derivadas(1)=-y-z;
derivadas(2)=x+0.2*y;
derivadas(3)=0.2+(x-5.7)*z;
end