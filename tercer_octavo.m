% ------------ parámetros refinados ------------
L  = 1;   c  = 1;
Ne = 20;  h  = L/Ne;
dt = 0.02;       % dt < h / c  (≈0.05)
T  = 2;   Nt = round(T/dt);

% ------------ matrices de elemento ------------
Me = h/6 * [2 1; 1 2];
Ke = c^2/h * [1 -1; -1 1];

% ------------ ensamblaje global ---------------
Nn  = Ne+1;
M = zeros(Nn);   K = zeros(Nn);
for e = 1:Ne
   idx = [e e+1];
   M(idx,idx) = M(idx,idx) + Me;
   K(idx,idx) = K(idx,idx) + Ke;
end

% ------------ aplicar Dirichlet ----------------
free = 2:Nn-1;
Mr = M(free,free);     Kr = K(free,free);

% versión lumped para abaratar la inversión
Mr = diag(sum(Mr,2));

A  = (dt^2) * (Mr\Kr);

% ------------ condiciones iniciales ------------
x  = linspace(0,L,Nn);
U0 = sin(pi*x(free)).';
U1 = U0 - 0.5*dt^2*(Mr\Kr*U0);   % arranque de orden 2

% ------------ integración explícita -----------
Unm1 = U0;  Un = U1;
Uhist = zeros(numel(free),Nt+1);
Uhist(:,1:2) = [U0,U1];

for n = 2:Nt
    Unp1 = 2*Un - Unm1 - A*Un;
    Uhist(:,n+1) = Unp1;
    Unm1 = Un;  Un = Unp1;
end

% ------------ reconstruir y graficar ----------
Ufull = [zeros(1,Nt+1); Uhist; zeros(1,Nt+1)];
tvec  = 0:dt:T;

figure
mesh(tvec,x,Ufull); view(45,30)
xlabel('t'), ylabel('x'), zlabel('u_h'), title('FEM + diferencias centradas')
