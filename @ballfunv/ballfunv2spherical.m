function varargout = ballfunv2spherical(v)
% Convert a ballfunv to spherical system and return Vr, Vlam, Vth

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% We might need to improve this code by doing n+2 and p+2 and interpolating the result

[Vx,Vy,Vz] = v.comp{:};
Vx = Vx.coeffs;
Vy = Vy.coeffs;
Vz = Vz.coeffs;

[m,n,p] = size(Vx);

% Useful spectral matrices
MsinL = trigspec.multmat(n, [0.5i;0;-0.5i] ); 
McosL = trigspec.multmat(n, [0.5;0;0.5] );
MsinT = trigspec.multmat(p, [0.5i;0;-0.5i] ); 
McosT = trigspec.multmat(p, [0.5;0;0.5] );

Vr = zeros(n,p,m);
Vlam = zeros(n,p,m);
Vth = zeros(n,p,m);

% Permute the scalars
Vx = permute(Vx,[2,3,1]);
Vy = permute(Vy,[2,3,1]);
Vz = permute(Vz,[2,3,1]);

% Loop over r
for k = 1:m
   Vr(:,:,k) = McosL*Vx(:,:,k)*MsinT.' + MsinL*Vy(:,:,k)*MsinT.' + Vz(:,:,k)*McosT.';
   Vlam(:,:,k) = -MsinL*Vx(:,:,k) + McosL*Vy(:,:,k);
   Vth(:,:,k) = McosL*Vx(:,:,k)*McosT.' + MsinL*Vy(:,:,k)*McosT.' - Vz(:,:,k)*MsinT.';
end

% Permute back
Vr = permute(Vr,[3,1,2]);
Vlam = permute(Vlam,[3,1,2]);
Vth = permute(Vth,[3,1,2]);

Vr = ballfun(Vr);
Vlam = ballfun(Vlam);
Vth = ballfun(Vth);

% Prepare output:
if ( nargout <= 1 ) 
    varargout = { [Vr,Vlam,Vth] };
else 
    varargout = { Vr,Vlam,Vth };
end
end