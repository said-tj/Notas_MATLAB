%======================================================================
% 1D Finite Element Method + Backward Euler para la ecuación de calor
%
%   du/dt = α * d²u/dx² + f(x,t),    0 ≤ x ≤ L,   0 ≤ t ≤ tf
%
% Condiciones de frontera:
%   u(0,t) = 0,    u(L,t) = 0
% Condición inicial:
%   u(x,0) = u0(x)
%
% Pasos:
%  1) Sustituir w(x,t) → aquí: aproximamos u(x,t) ≈ Σ Ni(x) U_i(t)
%  2) Ponderar (Galerkin) y obtener la forma débil
%  3) Integrar en cada elemento para formar M (masa) y K (rigidez)
%  4) Imponer condiciones iniciales y de frontera
%  5) Integrar en el tiempo (Backward Euler) y resolver el sistema lineal
%======================================================================

clear; clc;

%% Parámetros físicos y numéricos
L    = 1.0;      % Longitud del dominio en x
tf   = 0.5;      % Tiempo final
alpha = 1.0;     % Difusividad térmica

nx   = 10;       % número de elementos en x
dx   = L/nx;     % tamaño de elemento
x    = linspace(0,L,nx+1)';  % nodos: i=1..nx+1

nt   = 100;      % número de pasos temporales
dt   = tf/nt;    % tamaño de paso en t
t    = linspace(0,tf,nt+1);

%% Ensamblado de matrices de masa y rigidez
M = sparse(nx+1,nx+1);
K = sparse(nx+1,nx+1);

for e = 1:nx
    nodes = [e, e+1];
    h = dx;
    % Matrices locales (elemento lineal, peso Galerkin)
    Mloc = (h/6) * [2 1; 1 2];
    Kloc = (alpha/h) * [ 1 -1; -1 1];
    % Ensamblar al sistema global
    M(nodes,nodes) = M(nodes,nodes) + Mloc;
    K(nodes,nodes) = K(nodes,nodes) + Kloc;
end

%% Condiciones de frontera (Dirichlet homogéneas en extremos)
fixedDofs = [1, nx+1];
freeDofs  = setdiff(1:nx+1, fixedDofs);

%% Condición inicial
u0 = @(x) sin(pi*x);        % ejemplo: modo fundamental
U = u0(x);

%% Ensamblado del esquema en tiempo (Backward Euler)
A = M(freeDofs,freeDofs) + dt * K(freeDofs,freeDofs);

%% Fuente f(x,t) 
% (si no hay, déjala en cero; ejemplo genérico:)
f = @(x,t) 0.*x + 0.*t;   

%% Bucle temporal
for n = 1:nt
    tn = t(n+1);
    % Vector de carga en tiempo tn
    F_full = zeros(nx+1,1);
    for i = 1:nx+1
        F_full(i) = f(x(i), tn);
    end
    % Sólo liberamos los dofs libres
    b = M(freeDofs,freeDofs)*U(freeDofs) + dt * F_full(freeDofs);
    % Resolver sistema reducido
    U_new_free = A \ b;
    % Actualizar U, dejando en cero los nodos fijos
    U(freeDofs) = U_new_free;
    U(fixedDofs) = 0;
end

%% Resultados y gráfica
figure;
plot(x, U, '-o', 'LineWidth', 1.5);
xlabel('x'); ylabel('u(x,tf)');
title(sprintf('Solución aproximada en t = %.2f', tf));
grid on;
