%% Apply the 2/3 rule
function V = Orsag_Truncate(Vexp,S)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[~,n_tilde,p_tilde] = size(Vexp);
m = S(1); n = S(2); p = S(3);
[VexpX,VexpY,VexpZ] = Vexp.comp{:};
VexpX = VexpX.coeffs; VexpY = VexpY.coeffs; VexpZ = VexpZ.coeffs;
Vx = VexpX(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2),floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2));
Vy = VexpY(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2),floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2));
Vz = VexpZ(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2),floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2));
V = ballfunv(ballfun(Vx),ballfun(Vy),ballfun(Vz));
end