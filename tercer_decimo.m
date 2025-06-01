clear;  clc;  close all;

%% Primero. Parámetros básicos
L   = 1;          % dominio espacial [0,1]
c   = 1;          % velocidad de la onda
Ne  = 3;          % 3 sub-intervalos  → 4 nodos
h   = L/Ne;       % longitud de cada elemento
dt  = 0.2;        % paso temporal  (dt < h / c)
T   = 2.0;        % tiempo final
Nt  = round(T/dt);

%% Segundo. Matrices de elemento (P1)
Me = (h/6) * [2 1; 1 2];      % masa local
Ke = (1/h) * [1 -1; -1 1];    % rigidez local  (c² = 1)

%% Tercero. Ensamblaje global (4 × 4)
Nn = Ne + 1;                   % 4 nodos globales
M  = zeros(Nn);  K  = zeros(Nn);
for e = 1:Ne
    idx = [e e+1];             % nodos locales del elemento e
    M(idx,idx) = M(idx,idx) + Me;
    K(idx,idx) = K(idx,idx) + Ke;
end

%% Cuarto. Dirichlet homogéneo  →  se eliminan nodos 0 y L
free = 2:3;                    % nodos interiores (MATLAB 2 y 3)
Mr = M(free,free);
Kr = K(free,free);

%% Quinto. Matriz A del esquema explícito centrado
A = (dt^2) * (Mr\Kr);          %  A = dt²·Mr⁻¹·Kr

%% Sexto. Condiciones iniciales  (u = sin(πx),   u_t = 0)
U0 = [ sin(pi/3);  sin(2*pi/3) ];      % desplazamiento en t = 0
U1 = U0 - 0.5*A*U0;                    % paso t = dt (Taylor 2º orden)

Uh = zeros(numel(free), Nt+1);         % almacenará U en todos los pasos
Uh(:,1) = U0;        % t = 0
Uh(:,2) = U1;        % t = dt

%% Septimo. Bucle temporal  (diferencias centradas)
Unm1 = U0;  Un = U1;
for n = 2:Nt
    Unp1 = 2*Un - Unm1 - A*Un;   % esquema explícito
    Uh(:,n+1) = Unp1;
    Unm1 = Un;
    Un   = Unp1;
end

%% Octavo. Reconstrucción con nodos de contorno (u = 0) para graficar
tvec  = 0:dt:T;
Ufull = [ zeros(1,Nt+1);   % nodo x = 0
          Uh;              % nodos interiores
          zeros(1,Nt+1) ]; % nodo x = 1
xnod  = linspace(0,L,Nn);

%% Noveno. Solución analítica (para mostrar junto a la numérica)
[X,Tm] = meshgrid(xnod,tvec);
Uex = sin(pi*X) .* cos(pi*Tm);

%% Decimo Primero. Gráfica: superficie numérica + exacta
figure('Name','Onda 1-D: FEM (malla) y exacta (sombra)');
mesh(tvec, xnod, Ufull, 'FaceAlpha',1); hold on
mesh(tvec, xnod, Uex',  'EdgeColor','none','FaceAlpha',0.3);
view(45,25)
xlabel('t'); ylabel('x'); zlabel('u');
title('Solución Método Elementos Finitos (azul) y exacta (translúcida)');
grid on
