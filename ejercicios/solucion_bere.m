clear all; clc;

f_inferior = @(x,y) 2*x;
f_superior = @(x,y) 40;

a = 20; 
f_izquierda = @(x,y) 4*y;
f_derecha = @(x,y) 40 + 40 * sin(pi * y / 10);
b = 10; 
n = 5;
m = 6;
hx = a / m; 
hy = b / n;


A1 = diag(-2 * (1 + (hx/hy)^2) * ones(m-1,1)) + diag((hx/hy)^2 * ones(m-2,1),1) + diag((hx/hy)^2 * ones(m-2,1),-1);

A = zeros((m-1)*(n-1));

for i = 1:n-1
    A(((i-1)*(m-1)+1):(i*(m-1)), ((i-1)*(m-1)+1):(i*(m-1))) = A1;
end

for i = 1:n-2
    A(((i-1)*(m-1)+1):(i*(m-1)), (i*(m-1)+1):((i+1)*(m-1))) = eye(m-1);
    A((i*(m-1)+1):((i+1)*(m-1)), ((i-1)*(m-1)+1):(i*(m-1))) = eye(m-1);
end

b_vec = zeros((m-1)*(n-1),1);

for i = 1:m-1
    x_i = i * hx;
    b_vec(i) = f_inferior(x_i, 0); 
    b_vec((m-1)*(n-2) + i) = f_superior(x_i, b); 
end

for j = 1:n-1
    y_j = j * hy;
    b_vec(((j-1)*(m-1)+1)) = f_izquierda(0, y_j); 
    b_vec(j*(m-1)) = f_derecha(a, y_j); 
end
w = a / b;
% bordes
W = zeros(m+1, n+1);
W(:,:) = w; 

disp('Matriz A:');
disp(A);

disp('Vector b:');