function varargout = PTdecomposition(v)
% PTDECOMPOSITION Poloidal-toroidal decomposition of a BALLFUNV
%   PTDECOMPOSITION(v) is the poloidal-toroidal decomposition of the BALLFUNV v
%   The Poloidal-toroidal decomposition of a divergence-free BALLFUNV v 
%   is of the form V = P + T

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[Vx,Vy,Vz] = v.comp{:};

% Get the discretization
[m,n,p] = size(Vx);

% Spectral matrices
Mr = ultraS.multmat(m, [0;1], 0);
MsinL = trigspec.multmat(n, [0.5i;0;-0.5i]); 
McosL = trigspec.multmat(n, [0.5;0;0.5]);
DF1L = 1i*spdiags((-floor(n/2):floor(n/2))', 0, n,n);

MsinT = trigspec.multmat(p, [0.5i;0;-0.5i]);
Msin2T = trigspec.multmat(p, [-0.25;0;0.5;0;-0.25]);
McossinT = trigspec.multmat(p, [0.25i;0;0;0;-0.25i]);
Msin2cosT = trigspec.multmat(p, [-.125;0;.125;0;.125;0;-.125]);
Msin3T = trigspec.multmat(p, [-.125i;0;.375i;0;-.375i;0;.125i]);
DF1T = 1i*spdiags((-floor(p/2):floor(p/2))', 0, p, p);


% Extract coeffs
Vx = Vx.coeffs;
Vy = Vy.coeffs;
Vz = Vz.coeffs;

% Permute Vx, Vy and Vz
Vx = permute(Vx,[2,3,1]);
Vy = permute(Vy,[2,3,1]);
Vz = permute(Vz,[2,3,1]);

% Rhs of the poloidal part
RhsP = zeros(size(Vx));
% Rhs of the toroidal part
RhsT = zeros(size(Vx));

% Loop over r
for k = 1:m
   RhsP(:,:,k) = -McosL*Vx(:,:,k)*Msin3T.' - MsinL*Vy(:,:,k)*Msin3T.' - Vz(:,:,k)*Msin2cosT.';
   RhsT(:,:,k) = MsinL*Vx(:,:,k)*MsinT.'*DF1T.'*MsinT.' + DF1L*McosL*Vx(:,:,k)*McossinT.'...
               - McosL*Vy(:,:,k)*MsinT.'*DF1T.'*MsinT.' + DF1L*MsinL*Vy(:,:,k)*McossinT.'...
               - DF1L*Vz(:,:,k)*Msin2T.';
end

% Permute back
RhsP = permute(RhsP,[3,1,2]);
RhsT = permute(RhsT,[3,1,2]);

% Multiply RhsP by r
for k = 1:p
   RhsP(:,:,k) = Mr*RhsP(:,:,k);
end

% Poloidal and toroidal scalars
P = ballfunv.PTequation(ballfun(RhsP));
T = ballfunv.PTequation(ballfun(RhsT));

% Prepare output:
if ( nargout <= 1 ) 
    varargout = { [P,T] };
else 
    varargout = { P, T };
end
end
