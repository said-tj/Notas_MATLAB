% ------------------------------------------------------------
%  WAVE 1-D  (u_tt = c² u_xx)  ·  3 elementos lineales  ·  Δt = 0.2
%  Guarda la malla completa  Ugrid(nodo , tiempo)
% ------------------------------------------------------------
clc;  clear;

% ---- 1 ·  Parámetros básicos
L  = 1;            % dominio espacial [0,1]
c  = 1;            % velocidad
Ne = 3;            % nº de elementos  -> 4 nodos
h  = L/Ne;         % 1/3
dt = 0.2;          % paso temporal
T  = 2;            % tiempo final
Nt = round(T/dt);  % nº de pasos (0..Nt)
nodes = Ne+1;      % 4
x     = linspace(0,L,nodes)';   % coordenadas nodales
tvec  = 0:dt:T;                 % vector tiempo

% ---- 2 ·  Matrices de elemento (lineales)
Me = h/6 * [2 1; 1 2];
Ke = c^2/h * [ 1 -1 ; -1 1];

% ---- 3 ·  Ensamblaje global  (4×4)
M = zeros(nodes);  K = zeros(nodes);
for e = 1:Ne
    idx = [e e+1];
    M(idx,idx) = M(idx,idx) + Me;
    K(idx,idx) = K(idx,idx) + Ke;
end

% ---- 4 ·  Dirichlet homogéneo en nodos 1 y 4
free = 2:3;                 % DOF interiores
Mr  = M(free,free);
Kr  = K(free,free);
A   = dt^2 * (Mr\Kr);       % matriz dinámica   2×2

% ---- 5 ·  Condiciones iniciales
U0_int = sin(pi*x(free));   % u(x,0) = sin(pi x)   -> nodos 2 y 3
Um1_int = U0_int;           % u_t = 0  ->  U^{-1} = U^{0}

% ---- 6 ·  Almacenamiento completo de la malla
Ugrid = zeros(nodes, Nt+1);       % nodos × instantes
Ugrid(free,1) = U0_int;           % t = 0
% nodos 1 y 4 ya están en 0 (Dirichlet)

% ---- 7 ·  Integración en el tiempo
Un   = U0_int;
Unm1 = Um1_int;
for n = 1:Nt
    Unp1 = 2*Un - Unm1 - A*Un;   % esquema central
    % guardar instante n
    Ugrid(free,n+1) = Unp1;
    % avanzar
    Unm1 = Un;
    Un   = Unp1;
end

% ---- 8 ·  Resultados: matriz completa
disp('Matriz Ugrid  (filas: nodos 1..4  ·  columnas: t = 0,dt,...,T)');
disp(array2table(Ugrid,'VariableNames', ...
     compose('t=%.1f',tvec),'RowNames',compose('x=%.2f',x)));

% ---- 9 ·  Visualización rápida
figure
surf(tvec, x, Ugrid,'EdgeColor','k');
xlabel('t'), ylabel('x'), zlabel('u(x,t)')
title('Evolución  FEM (3 elem)  +  Δt = 0.2')
colormap jet, shading interp
