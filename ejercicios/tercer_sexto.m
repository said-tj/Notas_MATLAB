%% wave_heat_fem.m
%  FEM 1-D (espacio)  +  base cúbica en el tiempo
%  Construcción de la gran matriz M = kron(Exx , Lyy) ± …  que tu profesor pidió
clc; clear;

%% ▸ 1. Parámetros físicos y numéricos
L  = 30;     % longitud del dominio
tf = 2;      % tiempo final
n  = 1;      % nº de “elementos” espaciales     (malla muy grosera)
m  = 1;      % nº de pasos de tiempo            (idéntico comentario)
hx = L/n;    % paso espacial
ht = tf/m;   % paso temporal
alpha = 1;   % coeficiente que multiplica la 2ª derivada temporal

% *********************************************************************
% ▸ 2. MATRICES–BLOQUE 4×4 QUE HAY QUE TENER DEFINIDAS
%    (en tu código eran vectores random – aquí las ponemos coherentes)
%    · Exx  →   masa (o rigidez) en X
%    ·  E   →   masa identidad en X
%    · Lyy  →   rigidez (o masa) en Y      ——  lo decides según tu teoría
%    · Lt   →   masa/rigidez en T
% *********************************************************************

% --- Ejemplo muy simple con “hat” lineales en X y base cúbica en T -----
%     (Pon aquí los bloques exactos que te haya dado tu profesor.
%      Usamos algo genérico sólo para que el código funcione.)

% 4 nodos (n + 3 = 4)  →  base cúbica B-spline / C2
Exx = [  2 -2  0  0;
        -2  4 -2  0;
         0 -2  4 -2;
         0  0 -2  2 ] * hx/6;         % “masa” en X

E   = eye(4);                         % masa idéntica (si hace falta)

Lyy = [  6  -6   0   0;
        -6  12  -6   0;
         0  -6  12  -6;
         0   0  -6   6 ] / hx^3;      % “rigidez” en X (2ª deriv.)

Lt  = eye(4);                         % masa en el “tiempo” (sólo demo)

% ----------------------------------------------------------------------
% ▸ 3. ENSAMBLAJE GLOBAL MEDIANTE PRODUCTO DE KRONECKER
% ----------------------------------------------------------------------
dof_x = n + 3;       % = 4
dof_t = m + 3;       % = 4
Ntot  = dof_x * dof_t;

% Bloques Laplacianos / masa
M1 = kron( Exx , Lyy );     % parte espacial×espacial   (o Exx⊗Lyy)
M2 = kron( E   , Lt  );     % parte masa×masa           (o  E ⊗Lt)

% Matriz final; ajusta el signo según sea ecuación de Onda o de Calor
%   • Para la Onda          → M = M1 + alpha^2 * M2
%   • Para la Calor ( ∂t u) → M = M1 - alpha^2 * M2
M = M1 - alpha^2 * M2;

%% ▸ 4. Mostrar resultado
disp('Matriz global                  size =');
disp(size(M));
disp('----------------------------------------------');
disp('Sub-bloque 8×8 de la matriz final M :');
disp(M(1:8,1:8));           % solo una parte para no inundar la pantalla
