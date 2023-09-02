function xp = L1_SPL_1(x0, A, At, b, V, pdtol, pdmaxiter, cgtol, cgmaxiter)
b=b.*V;
for i = 1:size(A,2)
    A(:,i)=A(:,i).*V;
end
index=[];
for i=1:size(V,1)
    if V(i)==0
        index=[index,i];
    end
end
b(index,:)=[];
A(index,:)=[];
largescale = isa(A,'function_handle');
if (nargin < 6)
    pdtol = 1e-12;
end
if (nargin < 7)
    pdmaxiter = 1000;
end
if (nargin < 8)
    cgtol = 1e-12;
end
if (nargin < 9)
    cgmaxiter = 1000;
end
N = length(x0);
alpha = 0.01;
beta = 0.5;
mu = 10;
gradf0 = [zeros(N,1); ones(N,1)];
x = x0;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
fu1 = x - u;
fu2 = -x - u;
lamu1 = -1./fu1;
lamu2 = -1./fu2;
if (largescale)
  v = -A(lamu1-lamu2);
  Atv = At(v);
  rpri = A(x) - b;
else
  v = -A*(lamu1-lamu2);
  Atv = A'*v;
  rpri = (A*x - b);
end
sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*N/sdg;
rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
resnorm = norm([rdual; rcent; rpri]);
pditer = 0;
done = (sdg < pdtol) | (pditer >= pdmaxiter);
while (~done)
  pditer = pditer + 1; 
  w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;
  w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
  w3 = -rpri;
  sig1 = -lamu1./fu1 - lamu2./fu2;
  sig2 = lamu1./fu1 - lamu2./fu2;
  sigx = sig1 - sig2.^2./sig1;  
  if (largescale)
    w1p = w3 - A(w1./sigx - w2.*sig2./(sigx.*sig1));
    h11pfun = @(z) -A(1./sigx.*At(z));
    [dv, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('Primal-dual: Cannot solve system.  Returning previous iterate.');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - At(dv))./sigx;
    Adx = A(dx);
    Atdv = At(dv);
  else
    H11p = -A*diag(1./sigx)*A';
    w1p = w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1));
    [dv,hcond] = linsolve(H11p,w1p);
    if (hcond < 1e-128)
      disp('Primal-dual: Matrix ill-conditioned.  Returning previous iterate.');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - A'*dv)./sigx;
    Adx = A*dx;
    Atdv = A'*dv;
  end
  du = (w2 - sig2.*dx)./sig1;
  dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;
  dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2; 
  indp = find(dlamu1 < 0);  indn = find(dlamu2 < 0);
  s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);
  indp = find((dx-du) > 0);  indn = find((-dx-du) > 0);
  s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);  
  backiter = 0;
  xp = x + s*dx;  up = u + s*du; 
  vp = v + s*dv;  Atvp = Atv + s*Atdv; 
  lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
  fu1p = xp - up;  fu2p = -xp - up;  
  rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
  rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
  rpp = rpri + s*Adx;  
  if (isnan(fu1p(1))==1)   
    xp = x;
    return
  end
  while(norm([rdp; rcp; rpp]) > (1-alpha*s)*resnorm)
    s = beta*s;
    xp = x + s*dx;  up = u + s*du; 
    vp = v + s*dv;  Atvp = Atv + s*Atdv; 
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
    fu1p = xp - up;  fu2p = -xp - up;  
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    backiter = backiter+1;
    if (backiter > 32)
      disp('Stuck backtracking, returning last iterate.')
      xp = x;
      return
    end   
  end
  x = xp;  u = up;
  v = vp;  Atv = Atvp; 
  lamu1 = lamu1p;  lamu2 = lamu2p;
  fu1 = fu1p;  fu2 = fu2p;    
  sdg = -(fu1'*lamu1 + fu2'*lamu2);
  tau = mu*2*N/sdg;
  rpri = rpp;
  rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
  rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
  resnorm = norm([rdual; rcent; rpri]); 
  done = (sdg < pdtol) | (pditer >= pdmaxiter); 
  if (largescale)
    fprintf('CG Res = %8.3e, CG Iter = %d\n', cgres, cgiter);
  end  
end
