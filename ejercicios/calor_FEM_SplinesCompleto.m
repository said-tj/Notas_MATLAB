% =======================================================================
%   calor_FEM_SplinesCompleto.m
% -----------------------------------------------------------------------
%   Ecuación de calor 1-D:
%        u_t = α u_xx        0 < x < L , 0 < t ≤ T
%   Dirichlet: u(0,t)=g0(t) , u(L,t)=gL(t)
%   Inicial  : u(x,0)=f(x)
%
%   Discretización espacial:  SPLINES CÚBICOS BÁSICOS UNIFORMES
%   Discretización temporal: θ-método (θ=1 ⇒ Euler implícito,
%                                         θ=0.5 ⇒ Crank-Nicolson, …)
%
%   Imprime todo lo que exige la guía del profesor:
%     1) Sustitución y nº de grados de libertad
%     2) Ponderación Galerkin
%     3) Integrales  C_ij (=∫S_i S_j)  y  K_ij (=∫S_i' S_j')
%        + matrices M y K completas o recortadas
%     4) Condiciones iniciales / frontera
%     5) Sistema lineal (A, B, factor LU)  y evolución gráfica u(x,t)
%
%   Autor: Saïdt  ·  Junio-2025
% =======================================================================

clear;  clc;  close all;
format rational          % <- muestra fracciones cuando sea posible

%% ------------------------  DATOS DE ENTRADA  --------------------------
fprintf('==========  FEM 1-D con splines cúbicos  ==========\n');
L      = input('Longitud del dominio   L   = ');
alpha  = input('Coef. de difusión      α   = ');
N      = input('Número de *elementos*  N   = ');
dt     = input('Paso temporal          Δt  = ');
Tfin   = input('Tiempo final           T   = ');
theta  = input('θ del método (1=implícito, 0.5=C-N, etc.) = ');

% --- Funciones problema (modifícalas si hace falta) --------------------
f  = @(x) sin(pi*x/L);    % condición inicial u(x,0)
g0 = @(t) 0*t;            % Dirichlet en x = 0
gL = @(t) 0*t;            % Dirichlet en x = L

%% ------------------------  MALLA Y DOF  -------------------------------
h    = L/N;                       % paso
x    = linspace(0,L,N+1)';        % nodos  x_0 … x_N
nDOF = N-1;                       % nodos interiores

fprintf('\nPaso 1) Sustitución:  u ≈ Σ U_i(t) S_i(x)\n');
fprintf('        → nDOF = %d grados de libertad (sin fronteras)\n', nDOF);

fprintf('\nPaso 2) Ponderación:  v = S_k(x)  (Galerkin)\n');

%% ----------------  COEFICIENTES DEL PROFESOR  -------------------------
%   C_{|i-j|}  multiplican *h*          (matriz de masa)
c0 =  5/16;     c1 =  3/112;    c2 = 129/2240;   c3 = 1/112;

%   K_{|i-j|}  multiplican *1/h*        (matriz de rigidez)
k0 = 13/20;     k1 = -13/60;    k2 =   1/60;    k3 = 0;

%% ---------------  MATRICES DE MASA (M) Y RIGIDEZ (K)  -----------------
fprintf('\nPaso 3) Integrales C_ij y K_ij:\n');
fprintf('   |i-j|   C_ij (=∫S_i S_j)        K_ij (=∫S_i'' S_j'')\n');
fprintf('     0     %s·h          %s/h\n',rat(c0),rat(k0));
fprintf('     1     %s·h          %s/h\n',rat(c1),rat(k1));
fprintf('     2     %s·h          %s/h\n',rat(c2),rat(k2));
fprintf('     3     %s·h          %s/h\n',rat(c3),rat(k3));

% vectores para spdiags
mainC = c0*ones(nDOF,1);  off1C=c1*ones(nDOF,1);  off2C=c2*ones(nDOF,1);  off3C=c3*ones(nDOF,1);
mainK = k0*ones(nDOF,1);  off1K=k1*ones(nDOF,1);  off2K=k2*ones(nDOF,1);  off3K=k3*ones(nDOF,1);

M = h      * spdiags([off3C off2C off1C mainC off1C off2C off3C], -3:3, nDOF, nDOF);
K = alpha/h* spdiags([off3K off2K off1K mainK off1K off2K off3K], -3:3, nDOF, nDOF);

fprintf('\n   Matriz de masa M — tamaño %d×%d, banda 7\n',size(M));
if N <= 12, disp(full(M)); else, disp('   (N grande ⇒ omito impresión completa)'); end
fprintf('\n   Matriz de rigidez K — tamaño %d×%d, banda 7\n',size(K));
if N <= 12, disp(full(K)); else, disp('   (N grande ⇒ omito impresión completa)'); end

%% ----------------  CONDICIONES INICIALES / FRONTERA  ------------------
fprintf('\nPaso 4) Condiciones iniciales y Dirichlet (homogéneas aquí)\n');
U = f(x(2:end-1));            % vector inicial
disp('   U(t=0) ='); disp(U.');

%% ----------------  SISTEMA LINEAL EN EL TIEMPO  -----------------------
fprintf('\nPaso 5) θ-método:  (M+θΔtK)U^{n+1} = (M-(1-θ)ΔtK)U^{n}\n');
A = M + theta*dt*K;
B = M - (1-theta)*dt*K;
[Llu,Ulu] = lu(A);            % factor LU (una sola vez)
fprintf('   Factorización LU concluida.\n');

%% ----------------  BUCLE TEMPORAL + GRÁFICA  --------------------------
nSteps = ceil(Tfin/dt);   t=0;
figure('Name','u(x,t) - FEM splines cúbicos');
for n = 1:nSteps
    t   = t + dt;
    rhs = B*U;
    U   = Ulu \ (Llu \ rhs);
    % reconstruir con Dirichlet para graficar
    Uplot = [g0(t); U; gL(t)];
    plot(x,Uplot,'-o','LineWidth',1.4);
    xlabel('x'); ylabel('u(x,t)'); grid on;
    title(sprintf('Solución FEM      t = %.4f',t)); ylim([-1 1]);
    drawnow;
end
fprintf('\n>>> Simulación terminada sin errores.\n');
