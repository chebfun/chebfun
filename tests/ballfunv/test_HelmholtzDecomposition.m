function pass = test_HelmholtzDecomposition( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e6*pref.techPrefs.chebfuneps;

%% Two-component form

% Example 1 :
vx = ballfun(@(x,y,z)cos(x.*y));
vy = ballfun(@(x,y,z)sin(x.*z));
vz = ballfun(@(x,y,z)cos(y.*z));
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
w = grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi);
pass(1) = norm(v-w)<tol;

% Example 2 :
f = ballfun(@(x,y,z)cos(y.*z));
P = ballfun(@(x,y,z)cos(x.*y));
T = ballfun(@(x,y,z)sin(x.*z));
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
w = grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi);
pass(2) = norm(v-w)<tol;

%% Three-component form

% Example 3 :
vx = ballfun(@(x,y,z)cos(x.*y));
vy = ballfun(@(x,y,z)sin(x.*z));
vz = ballfun(@(x,y,z)cos(y.*z));
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
w = grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi); 
pass(3) = norm(v-w)<tol;

% Example 4 :
f = ballfun(@(x,y,z)cos(y.*z));
P = ballfun(@(x,y,z)cos(x.*y));
T = ballfun(@(x,y,z)sin(x.*z));
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
w = grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi); 
pass(4) = norm(v-w)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end