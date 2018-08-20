function pass = test_HelmholtzDecomposition()

%% Two-component form

S = [51,51,51];

% Example 1 :
vx = ballfun(@(x,y,z)cos(x.*y),'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)cos(y.*z),'cart',S);
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
pass(1) = isequal(v,grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi));

% Example 2 :
f = ballfun(@(x,y,z)cos(y.*z),'cart',S);
P = ballfun(@(x,y,z)cos(x.*y),'cart',S);
T = ballfun(@(x,y,z)sin(x.*z),'cart',S);
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
pass(2) = isequal(v,grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi));

S = [51,51,51];

% Example 3 :
vx = ballfun(@(x,y,z)cos(x.*y),'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)cos(y.*z),'cart',S);
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
pass(3) = isequal(v,grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi));

% Example 4 :
f = ballfun(@(x,y,z)cos(y.*z),'cart',S);
P = ballfun(@(x,y,z)cos(x.*y),'cart',S);
T = ballfun(@(x,y,z)sin(x.*z),'cart',S);
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
pass(4) = isequal(v,grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi));

S = [53,53,53];

% Example 5 :
vx = ballfun(@(x,y,z)cos(x.*y),'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)cos(y.*z),'cart',S);
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
pass(5) = isequal(v,grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi));

% Example 6 :
f = ballfun(@(x,y,z)cos(y.*z),'cart',S);
P = ballfun(@(x,y,z)cos(x.*y),'cart',S);
T = ballfun(@(x,y,z)sin(x.*z),'cart',S);
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi] = HelmholtzDecomposition(v);
pass(6) = isequal(v,grad(f)+ballfunv.PT2ballfunv(Ppsi,Tpsi));

%% Three-component form

S = [51,51,51];

% Example 7 :
vx = ballfun(@(x,y,z)cos(x.*y),'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)cos(y.*z),'cart',S);
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
pass(7) = isequal(v,grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi));

% Example 8 :
f = ballfun(@(x,y,z)cos(y.*z),'cart',S);
P = ballfun(@(x,y,z)cos(x.*y),'cart',S);
T = ballfun(@(x,y,z)sin(x.*z),'cart',S);
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
pass(8) = isequal(v,grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi));

S = [51,51,51];

% Example 9 :
vx = ballfun(@(x,y,z)cos(x.*y),'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)cos(y.*z),'cart',S);
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
pass(9) = isequal(v,grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi));

% Example 10 :
f = ballfun(@(x,y,z)cos(y.*z),'cart',S);
P = ballfun(@(x,y,z)cos(x.*y),'cart',S);
T = ballfun(@(x,y,z)sin(x.*z),'cart',S);
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
pass(10) = isequal(v,grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi));

S = [53,53,53];

% Example 11 :
vx = ballfun(@(x,y,z)cos(x.*y),'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)cos(y.*z),'cart',S);
v = ballfunv(vx,vy,vz);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
pass(11) = isequal(v,grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi));

% Example 12 :
f = ballfun(@(x,y,z)cos(y.*z),'cart',S);
P = ballfun(@(x,y,z)cos(x.*y),'cart',S);
T = ballfun(@(x,y,z)sin(x.*z),'cart',S);
v = grad(f) + ballfunv.PT2ballfunv(P,T);
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
pass(12) = isequal(v,grad(f)+curl(ballfunv.PT2ballfunv(Ppsi,Tpsi))+grad(phi));

if (nargout > 0)
    pass = all(pass(:));
end
end