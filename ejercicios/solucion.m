%% Programa MATLAB: Solución de la ecuación de Laplace en una placa por diferencias finitas

% Parámetros geométricos y de malla
x0 = 0;    % extremo izquierdo en x
xN = 20;   % extremo derecho en x
y0 = 0;    % extremo inferior en y
yM = 10;   % extremo superior en y

n = 5;  % número de subintervalos en x (de modo que se tienen n+1 puntos en x)
m = 5;  % número de subintervalos en y (m+1 puntos en y)

hx = (xN - x0)/n;    % paso en x, hx = 20/5 = 4
hy = (yM - y0)/m;    % paso en y, hy = 10/5 = 2

% Tolerancia y número máximo de iteraciones
TOL = 1e-4;
Nmax = 100;

% Constantes del método
lambda = (hx^2)/(hy^2);  % lambda = 16/4 = 4
mu = 2*(1+lambda);       % mu = 10

% Se generan los puntos interiores.
% NOTA: Los nodos de frontera se incorporan a los términos fuente de la iteración.
% Se definen los puntos interiores (excluyendo los nodos donde se aplican condiciones de frontera)
x = (hx : hx : xN - hx);   % puntos interiores en x: 4,8,12,16
y = (hy : hy : yM - hy);   % puntos interiores en y: 2,4,6,8

nx = length(x);  % número de puntos interiores en x (n-1 = 4)
ny = length(y);  % número de puntos interiores en y (m-1 = 4)

% Inicializar la matriz de solución en los nodos interiores
w = zeros(nx, ny);

% Definir las funciones de condiciones de frontera:
gIzq = @(y) 4*y;                      % u(0,y)
gDer = @(y) 40 + 40*sin(pi*y/10);       % u(20,y)
fInf = @(x) 2*x;                      % u(x,0)
fSup = @(x) 40;                       % u(x,10)

% Inicialización del contador de iteraciones y la norma
iter = 0;
norm_val = Inf;

%% Inicio de iteraciones Gauss-Seidel
while iter < Nmax
    iter = iter + 1;
    norm_val = 0;  % se reinicia la norma para cada barrido

    %%% Paso 7: Actualizar el nodo superior izquierdo interior (i = 1, j = ny)
    i = 1; j = ny;  
    % La fórmula de actualización (considerando f(x,y)=0):
    % z = [ gIzq(y_j) + lambda*fSup(x_i) + lambda*w(1,j-1) + w(2,j) ] / mu;
    z = ( gIzq( y(j) ) + lambda * fSup(x(i)) + lambda * w(i, j-1) + w(i+1, j) ) / mu;
    diff = abs(z - w(i,j));
    if diff > norm_val, norm_val = diff; end
    w(i,j) = z;

    %%% Paso 8: Actualizar los puntos en la fila superior interior para i = 2,...,nx-1 (con j = ny)
    for i = 2:(nx-1)
        z = ( lambda * fSup(x(i)) + w(i-1, ny) + w(i+1, ny) + lambda * w(i, ny-1) ) / mu;
        diff = abs(z - w(i,ny));
        if diff > norm_val, norm_val = diff; end
        w(i,ny) = z;
    end

    %%% Paso 9: Actualizar el nodo superior derecho interior (i = nx, j = ny)
    i = nx; j = ny;
    z = ( gDer( y(j) ) + lambda * fSup(x(i)) + w(i-1, ny) + lambda * w(i, ny-1) ) / mu;
    diff = abs(z - w(i,j));
    if diff > norm_val, norm_val = diff; end
    w(i,j) = z;

    %%% Paso 10: Para j = ny-1 hasta 2 (recorriendo hacia abajo)
    for j = (ny-1):-1:2
        % Paso 11: Actualizar el primer nodo de la fila (i = 1)
        i = 1;
        z = ( gIzq(y(j)) + lambda * w(i, j+1) + lambda * w(i, j-1) + w(i+1, j) ) / mu;
        diff = abs(z - w(i,j));
        if diff > norm_val, norm_val = diff; end
        w(i,j) = z;

        % Paso 12: Actualizar los nodos centrales de la fila para i = 2,...,nx-1
        for i = 2:(nx-1)
            z = ( w(i-1,j) + lambda * w(i, j+1) + w(i+1,j) + lambda * w(i, j-1) ) / mu;
            diff = abs(z - w(i,j));
            if diff > norm_val, norm_val = diff; end
            w(i,j) = z;
        end

        % Paso 13: Actualizar el último nodo de la fila (i = nx)
        i = nx;
        z = ( gDer(y(j)) + w(i-1,j) + lambda * w(i, j+1) + lambda * w(i, j-1) ) / mu;
        diff = abs(z - w(i,j));
        if diff > norm_val, norm_val = diff; end
        w(i,j) = z;
    end

    %%% Paso 14: Actualizar el nodo inferior izquierdo interior (i = 1, j = 1)
    % En este paso se incorporan las condiciones de frontera en x=0 y y=0:
    % z = [ gIzq(y_1) + lambda*fInf(x_1) + lambda*w(1,2) + w(2,1) ] / mu;
    i = 1; j = 1;
    z = ( gIzq( y(j) ) + lambda * fInf(x(i)) + lambda * w(i,2) + w(i+1, j) ) / mu;
    diff = abs(z - w(i,j));
    if diff > norm_val, norm_val = diff; end
    w(i,j) = z;

    %%% Paso 15: Actualizar los puntos de la fila inferior interior para i = 2,...,nx-1 (con j = 1)
    for i = 2:(nx-1)
        % Aquí se usa la condición inferior: fInf(x) = 2x.
        z = ( lambda * fInf(x(i)) + w(i-1,1) + lambda * w(i,2) + w(i+1,1) ) / mu;
        diff = abs(z - w(i,1));
        if diff > norm_val, norm_val = diff; end
        w(i,1) = z;
    end

    %%% Paso 16: Actualizar el nodo inferior derecho interior (i = nx, j = 1)
    i = nx; j = 1;
    z = ( gDer( y(j) ) + lambda * fInf(x(i)) + w(i-1,1) + lambda * w(i,2) ) / mu;
    diff = abs(z - w(i,j));
    if diff > norm_val, norm_val = diff; end
    w(i,j) = z;

    %%% Paso 17: Verificar la convergencia
    if norm_val <= TOL
        fprintf('Convergencia alcanzada en %d iteraciones.\n', iter);
        break;
    end
end

if iter >= Nmax && norm_val > TOL
    disp('Se excedió el número máximo de iteraciones.');
end

%% Mostrar resultados
% Los valores de w(i,j) aproximan u(x,y) en los nodos interiores (x_i, y_j)
disp('Resultados (nodos interiores):');
for i = 1:nx
    for j = 1:ny
        fprintf('u(%.2f, %.2f) ≈ %.6f\n', x(i), y(j), w(i,j));
    end
end



%% Supongamos que ya tienes definida la solución interior 'w' de dimensión (n-1) x (m-1)
% y además los parámetros y funciones de frontera:
% Definición del dominio
x0 = 0;    xN = 20;
y0 = 0;    yM = 10;

n = 5;  % número de subintervalos en x  =>  n+1 = 6 puntos en x
m = 5;  % número de subintervalos en y  =>  m+1 = 6 puntos en y

hx = (xN - x0) / n;    % hx = 4
hy = (yM - y0) / m;    % hy = 2

% Funciones de condiciones de frontera:
gIzq = @(y) 4*y;                      % u(0,y)
gDer = @(y) 40 + 40*sin(pi*y/10);       % u(20,y)
fInf = @(x) 2*x;                      % u(x,0)
fSup = @(x) 40;                       % u(x,10)

%% Construir la malla completa (incluyendo fronteras)
% Crear vectores con todos los puntos en x e y
x_full = linspace(x0, xN, n+1);  % 0,4,8,12,16,20
y_full = linspace(y0, yM, m+1);  % 0,2,4,6,8,10

% Inicializar la matriz de solución completa U
U = zeros(m+1, n+1);  % Usamos (m+1) filas y (n+1) columnas

%% Asignar los valores en la frontera
% Borde inferior: y = y0
for i = 1:length(x_full)
    U(1,i) = fInf(x_full(i));
end

% Borde superior: y = yM
for i = 1:length(x_full)
    U(end,i) = fSup(x_full(i));
end

% Borde izquierdo: x = x0
for j = 1:length(y_full)
    U(j,1) = gIzq(y_full(j));
end

% Borde derecho: x = xN
for j = 1:length(y_full)
    U(j,end) = gDer(y_full(j));
end

%% Incorporar los valores del interior obtenidos (matriz w)
% Los valores interiores de w corresponden a los nodos: 
% filas: j = 2 : m   y columnas: i = 2 : n
% (recordar que en MATLAB la primera fila corresponde a y=y0 y la última a y=yM)
[nInt, mInt] = size(w);  % se espera que nInt = n-1 y mInt = m-1

for i = 1:nInt
    for j = 1:mInt
        % Nota:
        %   - La variable i de w corresponde al nodo en x en la posición i+1 en U.
        %   - La variable j de w corresponde al nodo en y en la posición j+1 en U.
        U(j+1, i+1) = w(i,j);
    end
end

%% Graficar la distribución de temperatura
% Generar la malla en 2D para el gráfico
[X, Y] = meshgrid(x_full, y_full);

figure;
surf(X, Y, U);         % Gráfico 3D de la superficie
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Distribución de temperatura en la placa');
shading interp;        % Interpolación de la superficie para mejor visualización
colorbar;              % Muestra la escala de colores

% También puedes graficar un contorno para ver las líneas de nivel
figure;
contourf(X, Y, U, 20); % 20 niveles de contorno
xlabel('x');
ylabel('y');
title('Contorno de la distribución de temperatura');
colorbar;
