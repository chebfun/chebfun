function varargout = HelmholtzDecomposition(v)
%HELMHOLTZDECOMPOSITION Helmholtz decomposition of a BALLFUNV
%   [F,P,T] = HELMHOLTZDECOMPOSITION(V) returns the two-component form of
%   the Helmholtz decomposition, where 
%   V = grad(F) + curl(curl(rP)) + curl(rT).
%   r denotes the vector of length sqrt(x^2+y^2+z^2) in the radial
%   direction.
%
%   [F,P,T,Phi] = HELMHOLTZDECOMPOSITION(V) returns the three-component form of
%   the Helmholtz decomposition, where
%   V = grad(F) + curl(curl(rP)) + curl(rT) + grad(phi).
%
%   See also PTdecomposition.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check if one component is empty
if isempty(v)
    isVempty = 1;
else
    isVempty = cellfun(@isempty, v.comp, 'UniformOutput', false);
    isVempty = any(cell2mat(isVempty));
end
if isVempty
    error('BALLFUNV:HelmholtzDecomposition:input', ...
          'ballfunv must not have an empty component');    
end

if nargout == 3
    % Compute the decomposition
    [f,Ppsi,Tpsi] = HelmholtzDecomposition_2(v);
    % Return the decomposition
    varargout = { f,Ppsi,Tpsi };
    
elseif nargout == 4
    % Compute the decomposition
    [f,Ppsi,Tpsi,phi] = HelmholtzDecomposition_3(v);
    % Return the decomposition
    varargout = { f,Ppsi,Tpsi,phi };
    
else
    error('BALLFUNV:HelmholtzDecomposition:unknown', ...
          'Undefined function ''HelhmoltzDecomposition'' for %d output arguments', nargout);    
end
end

function varargout = HelmholtzDecomposition_2(v)
% Compute the two-component form of the Helmholtz decomposition

% Compute the divergence of v
div_v = div(v);

% Size of div(v)
[m,n,p] = size(div_v);

% Increase the size if it's too small
m = max(m,5); n = max(n,5); p = max(p,5);

% Compute the boundary of the vector field v
v_Boundary = ComputeNormalBoundary(v,n,p);

% Solve the Poisson equation Delta f = div(v) with Boundary conditions
% df/dr = vr[r=1
f = helmholtz(div_v,0,v_Boundary,m,n,p,'neumann');

% Divergence-free vector field
v_1 = v - grad(f);

% PT decomposition of the vorticity
[Pv1,Tv1] = PTdecomposition(v_1);

% Return the decomposition
% Prepare output:
varargout = { f,Pv1,Tv1 };
end

function varargout = HelmholtzDecomposition_3(v)
% Compute the three-component form of the Helmholtz decomposition

% Get the discretization
div_v = div(v);
[m,n,p] = size(div_v);

% Increase the size if it's too small
m = max(m,5); n = max(n,5); p = max(p,5);

% Solve the Poisson equation Delta f = div(v)
f = helmholtz(div_v,0,zeros(n,p),m,n,p);

% Divergence free vector field
v_1 = v - grad(f);

% Get the discretization
S = max(size(v_1),[],1);
m = S(1); n = S(2); p = S(3);

% Increase the size if it's too small
m = max(m,5); n = max(n,5); p = max(p,5);

% Compute the boundary of the vector field v_1
v_Boundary = ComputeNormalBoundary(v_1,n,p);

% Solve Delta phi = 0 with Neumann boundary condition
% r.grad(phi) = dphi/dr = r.v_1
zero = ballfun(0);
phi = helmholtz(zero,0,v_Boundary,m,n,p,'neumann');

% Divergence-free and tangential-free vector field
v_2 = v_1 - grad(phi);

% PT decomposition of the vorticity
[Pv,Tv] = PTdecomposition(v_2);

% Use the PT decomposition of v to compute phi
% Solve 2 Poisson equations to find Ppsi and Tpsi
[m,n,p] = size(Tv);

% Increase the size if it's too small
m = max(m,5); n = max(n,5); p = max(p,5);

Ppsi = helmholtz(-Tv,0,zeros(n,p),m,n,p);
Tpsi = Pv;

% Return the decomposition
varargout = { f,Ppsi,Tpsi,phi };
end

function v_Boundary = ComputeNormalBoundary(v,n,p)
% COMPUTENORMALBOUNDARY compute vr|r=1

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the discretization
S = size(v);
mx = S(1,1); my = S(2,1); mz = S(3,1);

% Get the components
V = v.comp;
Vx = coeffs3(V{1},mx,n,p);
Vy = coeffs3(V{2},my,n,p);
Vz = coeffs3(V{3},mz,n,p);

% Evaluate at the boundary r = 1
Vx = reshape(sum(Vx,1),n,p);
Vy = reshape(sum(Vy,1),n,p);
Vz = reshape(sum(Vz,1),n,p);

% Useful spectral matrices
MsinL = trigspec.multmat(n, [0.5i;0;-0.5i] ); 
McosL = trigspec.multmat(n, [0.5;0;0.5] );
MsinT = trigspec.multmat(p, [0.5i;0;-0.5i] ); 
McosT = trigspec.multmat(p, [0.5;0;0.5] );

% Compute the boundary Vr|r=1
v_Boundary = McosL*Vx*MsinT.' + MsinL*Vy*MsinT.' + Vz*McosT.';
end