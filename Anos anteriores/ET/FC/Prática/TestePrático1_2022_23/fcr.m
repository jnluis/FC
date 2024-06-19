function F = fcr(xv,xold,yold,zold,h)
    F(1) = xv(1)-xold-(-yold-zold)*h;
    F(2) = xv(2)-yold-(xold-0.2*yold)*h;
    F(3) = xv(3)-zold-(0.2+(xold-1.8)*zold)*h;
end