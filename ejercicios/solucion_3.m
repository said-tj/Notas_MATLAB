a = 20;
b = 10;
f1 = @(x) 2*x;
f2 = @(x) 40 + 0 * x;
g1 = @(y) 4* y;
g2 = @(y) 40 + 40 * sin(pi*y / 10);
n = 5;
m = 6;
hx = a / n;
hy = b / m;
azul = (hy / hx)^z;

As = diag(azul * ones (m-1,1));
Ad = diag(ones (m-2))