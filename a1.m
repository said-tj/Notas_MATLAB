%% FEM – Heat equation – Paso 1: Sustituir w(x,t)
% Requiere Symbolic Math Toolbox

syms x t h alpha positive                     % símbolos globales
% Temperaturas nodales como funciones de t
syms T1(t) T2(t)

% Funciones de forma lineales en el elemento [0,h]
N1 = 1 - x/h;
N2 = x/h;

% Solución aproximada en el elemento
u_h = N1*T1 + N2*T2;                          

% Residuo de la ecuación de calor: R = u_t - alpha*u_xx
R = diff(u_h, t) - alpha*diff(u_h, x, 2);

% ===========================
% Paso 1: residuo ponderado
% ===========================
% Integramos w_i * R sobre el dominio del elemento
R1 = int(N1 * R, x, 0, h);   % Para w = N1
R2 = int(N2 * R, x, 0, h);   % Para w = N2

% Simplificamos y mostramos
R1 = simplify(R1);
R2 = simplify(R2);

disp('Residuo para el nodo 1 (N1):'); pretty(R1)
disp('Residuo para el nodo 2 (N2):'); pretty(R2)
