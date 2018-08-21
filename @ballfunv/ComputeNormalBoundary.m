function v_Boundary = ComputeNormalBoundary(v)
% COMPUTENORMALBOUNDARY compute vr|r=1

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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