function varargout = HelmholtzDecomposition(v)
% HELMHOLTZDECOMPOSITION Helmholtz decomposition of a BALLFUNV
%   HELMHOLTZDECOMPOSITION returns the two-component form or the three-component
%   form according to the number of output arguments
%   two-component form:
%   v = grad(f) + Psi, where Psi = curl(curl(rPpsi)) + curl(rTpsi)
%   three-component form:
%   v = grad(f) + curl(Psi) + grad(phi), where Psi = curl(curl(rPpsi)) + curl(rTpsi)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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

% Compute the boundary of the vector field v
v_Boundary = ComputeNormalBoundary(v);

% Solve the Poisson equation Delta f = div(v) with Boundary conditions
% df/dr = vr[r=1
div_v = div(v);
f = helmholtz_neumann(div_v,0,v_Boundary);


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
[m,n,p] = size(v);

% Solve the Poisson equation Delta f = div(v)
f = helmholtz(div(v),0,zeros(n,p));

% Divergence free vector field
v_1 = v - grad(f);

% Compute the boundary of the vector field v_1
v_Boundary = ComputeNormalBoundary(v_1);

% Solve Delta phi = 0 with Neumann boundary condition
% r.grad(phi) = dphi/dr = r.v_1
phi = helmholtz_neumann(ballfun(zeros(m,n,p)),0,v_Boundary);

% Divergence-free and tangential-free vector field
v_2 = v_1 - grad(phi);

% PT decomposition of the vorticity
[Pv,Tv] = PTdecomposition(v_2);

% Use the PT decomposition of v to compute phi
% Solve 2 Poisson equations to find Ppsi and Tpsi
Ppsi = helmholtz(-Tv,0,zeros(n,p));
Tpsi = Pv;

% Return the decomposition
varargout = { f,Ppsi,Tpsi,phi };
end

function v_Boundary = ComputeNormalBoundary(v)
% Compute Vr|r=1

% Get the discretization
[~,n,p] = size(v);

% Get the components
[Vx,Vy,Vz] = v.comp{:};
Vx = Vx.coeffs;
Vy = Vy.coeffs;
Vz = Vz.coeffs;

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