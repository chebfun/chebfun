function C = minus(A, B)
%-   CHEBOP minus.
%   A - B subtracts the operators of the CHEBOPs A and B
%
%   C = MINUS(A, B) is called for the syntax 'A - B'.
%
% Example:
%
%   L = chebop(@(u) -diff(u,2), [0,pi], 'dirichlet');
%   I = chebop(@(u) u, [0,pi]);
%   eigs(L-I)
%
% See CHEBOP/PLUS, CHEBOP/MTIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


C = plus(A, uminus(B), '');
    
end
