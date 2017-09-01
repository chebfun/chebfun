function H = cross( F, G )
%CROSS   Vector cross product.
%   CROSS(F, G) returns the cross product of the DISKFUNV objects F and G. 
%   If F and G both have two components, then it returns the DISKFUNV 
%   representing
%       CROSS(F,G) = F(1) * G(2) - F(2) * G(1)
%   where F = (F(1); F(2)) and G = (G(1); G(2)). 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    H = diskfunv();
    return
end

% Get components: 
Fc = F.components; 
Gc = G.components; 

% Do curl: 
H = Fc{1} .* Gc{2} - Fc{2} .* Gc{1};   % 2D cross
    
end


