% ---------------------------------------------------------------
%  B‐splines cúbicos  S_i(x)  i=0..5   (n = 3  -> n+3 = 6)
%  B‐splines cúbicos  S_j(t)  j=0..3   (m = 1  -> m+3 = 4)
%  Armado del sistema  (Kx ⊗ Mt  -  Mx ⊗ Kt) C = 0
%  y eliminación de las  filas/columnas que resultan nulas.
% ---------------------------------------------------------------
clc; clear;

%% 1 · Parámetros de malla
L  = 30;          % longitud espacial
tf = 2;           % tiempo final
n  = 3; m = 1;    % n = 3 → 6 S_i     ,  m = 1 → 4 S_j
hx = L / n;       % 10
ht = tf / m;      % 2
alpha = 1;        % c^2  (en tu notación)

%% 2 · Constantes de superposición  |i-r| <= 3
mu  = [ 4/5 , 3/10 , 1/20 , 1/120 ];              %  h  factor afuera
kap = [ 12/5 , -6/5 , 1/10 , -1/60 ];             %  1/h factor afuera

%% 3 · Construye Mx , Kx   (6 x 6)
Nx = n+3;  Nt = m+3;                 % 6 y 4
Mx = zeros(Nx);  Kx = zeros(Nx);

for i = 1:Nx
  for r = max(1,i-3):min(Nx,i+3)
        k = abs(i-r);
        Mx(i,r) = hx * mu(k+1);          % mu_0 .. mu_3
        Kx(i,r) = (1/hx) * kap(k+1);     % kap_0 .. kap_3
  end
end

%% 4 · Construye Mt , Kt   (4 x 4)
Mt = zeros(Nt);  Kt = zeros(Nt);

for j = 1:Nt
  for s = max(1,j-3):min(Nt,j+3)
        k = abs(j-s);
        Mt(j,s) = ht * mu(k+1);
        Kt(j,s) = (1/ht) * kap(k+1);
  end
end

%% 5 · Sistema completo 24 x 24
A = kron(Kx, Mt) - alpha * kron(Mx, Kt);    % 24 × 24

fprintf('Matriz original A:  tamaño %d × %d,  nnz = %d\n', ...
         size(A,1), size(A,2), nnz(A));

%% 6 · Detecta y elimina filas/columnas nulas
rows2keep = find(sum(abs(A),2)~=0);
Ared = A(rows2keep, rows2keep);

fprintf('Matriz reducida  :  tamaño %d × %d,  nnz = %d\n', ...
         size(Ared,1), size(Ared,2), nnz(Ared));

%% 7 · Visualiza la estructura esparsa
figure, spy(A),    title('Estructura de A (24 × 24)')
figure, spy(Ared), title('Estructura tras poda (16 × 16)')

%% 8 · (Opcional) Resuelve Ared * c = 0  con CI / BC
%  -> En este punto añades las ecuaciones de
%     • Dirichlet en x = 0,L  (colapsan los C_{0j} y C_{5j})
%     • Cond. inicial u(x,0)=f(x)  (fija la fila j = 0)
%     • Cond. inicial u_t(x,0)=g(x) (mezcla j = 0,1)
%     Aquí sólo mostramos la construcción del núcleo homogéneo.
