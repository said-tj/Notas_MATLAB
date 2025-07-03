clc, clear all;
n = 6;
m = 200;
g = @(x) 80*sin(pi * x/30);
L = 30;
hx = L/n;
tf = 200;
ht = tf/m;
x = linspace(0,L,n+1);
t = linspace(0,tf,m+1);
w = zeros(n+1,m+1);
for i = 2:n % Condición inicial
    w(i,1) = g(x(i));
end
for j = 2:m
    for i = 2:n
        w(i,j+1) = ht/(hx)^2 * (w(i-1,j)-2*w(i,j)+w(i+1,j)) + w(i,j);
    end
end
disp(w)