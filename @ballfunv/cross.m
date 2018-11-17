function H = cross(F,G)
%CROSS  Cross product of two BALLFUNVs.
%   CROSS(F, G) is the cross product of the BALLFUNV F and G.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the components
Fc = F.comp;
Gc = G.comp;

% Do cross: 
H = [ Fc{2} .* Gc{3} - Fc{3} .* Gc{2} ; ...
      Fc{3} .* Gc{1} - Fc{1} .* Gc{3} ; ...
      Fc{1} .* Gc{2} - Fc{2} .* Gc{1} ];
end