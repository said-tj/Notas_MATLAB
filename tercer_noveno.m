% ------------------------------------------------------------
%  Onda 1-D  :   u_tt = c^2 u_xx     (con c = 1 y L = 1)
%  Método de Elementos Finitos (3 elementos lineales, P1)
%  + Diferencias finitas centradas en el tiempo (explícito, 2.º orden)
%
%  Este script compara la solución numérica FEM (malla) con la
%  solución analítica   u(x,t) = sin(pi x) cos(pi t)
% ------------------------------------------------------------

clear;   clc;   close all;

%% ───────────────────── 1.  Parámetros del problema ─────────────────────
L   = 1;          % longitud física del dominio  x ∈ [0,L]
c   = 1;          % velocidad de la onda (aquí vale 1 => c^2 = 1)

Ne  = 3;          % número de sub-intervalos en X  (→ 4 nodos)
h   = L / Ne;     % tamaño de cada elemento  h = 1/3

dt  = 0.2;        % paso de tiempo  (cumple CFL: dt < h / c ≈ 0.333)
T   = 2.0;        % tiempo total de simulación
Nt  = round(T / dt);   % número de pasos de tiempo enteros

%% ───────────────────── 2.  Matrices de elemento (lineales) ─────────────
%  Para elementos P1:  masa (Me) y rigidez (Ke) de un elemento de longitud h
Me = (h/6) * [ 2  1 ;         % masa consistente 2×2
               1  2 ];

Ke = (1/h) * [ 1 -1 ;         % rigidez local 2×2  (c^2 = 1)
              -1  1 ];

%% ───────────────────── 3.  Ensamblaje global (4×4) ─────────────────────
Nn = Ne + 1;                   % número total de nodos  (4)
M  = zeros(Nn);                % matriz de masa global
K  = zeros(Nn);                % matriz de rigidez global

for e = 1:Ne
    idx = [e e+1];             % nodos locales del elemento e
    % sumar las contribuciones de cada elemento a las sub-matrices globales
    M(idx,idx) = M(idx,idx) + Me;
    K(idx,idx) = K(idx,idx) + Ke;
end

%% ───────────────────── 4.  Aplicar Dirichlet homogéneo ────────────────
%  Nodos de contorno:  x = 0 (nodo 1) y x = L (nodo 4) → u = 0
free = 2:3;                    % nodos interiores que quedan libres
Mr = M(free,free);             % masa reducida  (2×2)
Kr = K(free,free);             % rigidez reducida (2×2)

%% ───────────────────── 5.  Matriz A del esquema explícito ─────────────
%  U^{n+1} = 2U^n − U^{n−1} − dt^2 * (M_r^{-1} K_r) * U^n
A = (dt^2) * (Mr \ Kr);        %  dt^2 * Mr^{-1} Kr

%% ───────────────────── 6.  Condiciones iniciales ───────────────────────
%  u(x,0)  =  sin(pi x)
%  u_t(x,0) =  0
U0 = [ sin(pi/3) ; sin(2*pi/3) ];        % valores en nodos interiores

% Paso t = dt usando desarrollo de Taylor de orden 2 para mantener la
% precisión temporal (velocidad inicial cero ⇒ término lineal se anula)
U1 = U0 - 0.5 * A * U0;        % U^1 = U^0 − ½ dt^2 M^{-1}K U^0

% Almacén de resultados:  Uh(:,k) guarda solución en t = (k-1)*dt
Uh = zeros(numel(free), Nt + 1);
Uh(:,1) = U0;                  % t = 0
Uh(:,2) = U1;                  % t = dt

%% ───────────────────── 7.  Bucle temporal (diferencias centradas) ─────
Unm1 = U0;     % U^{n−1}
Un   = U1;     % U^{n}
for n = 2:Nt
    Unp1 = 2*Un - Unm1 - A*Un;    % esquema explícito centrado
    Uh(:,n+1) = Unp1;             % guardar resultado
    % avanzar los buffers
    Unm1 = Un;
    Un   = Unp1;
end

%% ───────────────────── 8.  Reconstrucción completa (incl. contorno) ───
%  Añadimos ceros en los nodos de contorno para trazar la superficie
tvec  = 0:dt:T;                         % vector de tiempos
Ufull = [ zeros(1,Nt+1) ;               % nodo x = 0
          Uh           ;               % nodos interiores
          zeros(1,Nt+1) ];             % nodo x = L
xnod  = linspace(0, L, Nn);             % nodos espaciales (0,1/3,2/3,1)

%% ───────────────────── 9.  Solución analítica (para referencia) ───────
%  u(x,t) = sin(pi x) * cos(pi t)
[X,Tm] = meshgrid(xnod, tvec);          % malla de nodos (t,x)
Uex = sin(pi*X) .* cos(pi*Tm);

%% ───────────────────── 10.  Gráfica de resultados ─────────────────────
figure('Name','Onda 1-D: FEM vs exacta');

% (a) Superficie numerica y exacta
subplot(1,2,1)
mesh(tvec, xnod, Ufull, 'FaceAlpha',1);  hold on
mesh(tvec, xnod, Uex',  'EdgeColor','none','FaceAlpha',0.3);
view(45,25)
xlabel('t');  ylabel('x');  zlabel('u');
title('Superficie FEM (malla)  y  exacta (sombra)');

% (b) Error L2 en función del tiempo
err = sqrt( h * sum( (Ufull - Uex').^2 ,1 ) );   % ∫|e|^2 dx ≈ h Σ nodos
subplot(1,2,2)
plot(tvec, err, 'k-o', 'LineWidth', 1.4);
xlabel('t');  ylabel('||error||_{L2}');
title('Error L_2  vs.  tiempo');
grid on;
