%% wave_fem_profesor.m
%  Modelo de prueba: ecuación de onda 1-D, Dirichlet homogéneo
%  malla: n+3 nodos (=4 si n=1)  –––   base cúbica en el tiempo (dof_t = m+3)

clc; clear; close all;

%% ───────────────────────── 1. Parámetros físico-numéricos ───────────────────────
L     = 30;          % longitud física del dominio
tf    =  2;          % tiempo final
c     =  1;          % velocidad de la onda  (cambia a alpha si es calor)

n     = 1;           % nº de “elementos” en el espacio  (muy tosco: n+3=4 nodos)
m     = 40;          % nº de pasos de tiempo (elige m suficientemente grande)
hx    =  L / n;      % paso espacial    (≃ 30   si n=1)
ht    =  tf / m;     % paso temporal
alpha = 1;           % aparece como α² en M = M1 – α²M2

dof_x = n + 3;       %  nodos en X (4)
dof_t = m + 3;       %  nodos en T (m+3)

%% ───────────────────────── 2. Bloques 4×4 de ensayo ────────────────────────────
%  ► Sustituye estas matrices por las que tu profesor haya entregado ◄
Exx = (hx/6)*[ 2 -2  0  0;
              -2  4 -2  0;
               0 -2  4 -2;
               0  0 -2  2 ];          % masa en X
E   = eye(dof_x);                     % identidad de apoyo
Lyy = (1/hx)*[ 1 -1  0  0;
              -1  2 -1  0;
               0 -1  2 -1;
               0  0 -1  1 ];          % rigidez en X
Lt  = eye(dof_t);                     % masa en T (para ejemplo)

%% ───────────────────────── 3. Matriz global M = M1 – α²M2 ───────────────────────
M1 = kron(Exx , Lyy);     % parte espacial × espacial
M2 = kron( E  , Lt );     % parte masa × masa
M  = M1 - alpha^2 * M2;

fprintf('\nMatriz global M (%d×%d) ensamblada.\n',size(M,1),size(M,2));

%% ───────────────────────── 4. Malla física y condiciones iniciales ─────────────
x_nodes = linspace(0,L,dof_x);        % nodos espaciales   (0, 30)
t_nodes = linspace(0,tf,m+1);         % tiempos que vamos a almacenar

% Grados de libertad interiores en X ⇒ i = 2 : dof_x-1  (porque i=1 y i=dof_x son Dirichlet 0)
intDOF = 2:(dof_x-1);                 % aquí sólo 2 y 3 cuando n=1

U0 = sin(pi*x_nodes(intDOF)/L).';     % u(x,0)  en nodos interiores
Ut0 = zeros(size(U0));                % velocidad inicial   du/dt(x,0) = 0
dof_int = numel(intDOF);              % =2 en este ejemplo

% ––– matriz de masa “lumped” de tamaño dof_int × dof_int para avance explícito
M_lumped = (hx/2)*diag([1 1]);        % (para 2 DOF interiores)
Minv     = diag(1./diag(M_lumped));   % inversa rápida

%% ───────────────────────── 5. Paso de tiempo: diferencias centradas ────────────
% Se usa el mismo esquema que en la respuesta anterior pero con L=30
dt = ht;                              % comodidad

% Condición CFL prudente  (h/c ≃ 30  → dt pequeño si n=1)
if dt > hx/c
     warning('dt = %.4g mayor que h/c = %.4g ; el esquema podría diverger.\n',...
              dt, hx/c);
end

% Primer paso (Taylor)  U1 = U0 – ½ dt² M⁻¹ K U0   (K = Lyy reducido)
K_red  = (1/hx)*[ 1 -1; -1 1 ];       % rigidez reducida (2×2)

U_hist      = zeros(dof_int , m+1);   % almacena solución en cada t
U_hist(:,1) = U0;
U1          = U0 - 0.5*dt^2*(Minv*K_red*U0);   % Ut0 = 0  →  término lineal ✗
U_hist(:,2) = U1;

% Matriz A de avance explícito:  U^{n+1} = (2I – dt² M⁻¹K) Uⁿ – U^{n-1}
A = 2*eye(dof_int) - dt^2 * (Minv*K_red);

for nStep = 2:m
    U_next = A*U_hist(:,nStep) - U_hist(:,nStep-1);
    U_hist(:,nStep+1) = U_next;
end

%% ───────────────────────── 6. Presentación del resultado ───────────────────────
fprintf('\nVector solución en el instante final t = %.2f s:\n',t_nodes(end));
disp(U_hist(:,end));

% Construimos la solución completa (añadiendo ceros en nodos de contorno)
U_full = [ zeros(1,m+1) ;             % x = 0
           U_hist ;                   % DOF interiores
           zeros(1,m+1) ];            % x = L

% Solución analítica para comparar
[Xgrid,Tgrid] = meshgrid(x_nodes , t_nodes);
U_exact = sin(pi*Xgrid'/L) .* cos(pi*c*Tgrid'/L);

%% ───────────────────────── 7.  Gráficas ─────────────────────────────────────────
figure('Name','Propagación de la onda FEM vs exacta','NumberTitle','off');

% (a) Superficie numérica
subplot(1,2,1);
surf(Tgrid , Xgrid , U_full,'EdgeColor','none');
xlabel('t'); ylabel('x'); zlabel('u_h');
title('Aproximación FEM (superficie)');
view(45,25);

% (b) Error L2 en el tiempo
errL2 = sqrt( hx * sum( (U_full - U_exact).^2 ,1 ) );
subplot(1,2,2);
plot(t_nodes , errL2 ,'k-o','LineWidth',1.4);
xlabel('t'); ylabel('\|error\|_{L2}');
title('Error L_2 frente a la solución exacta');
grid on;

