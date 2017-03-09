function C = minus(A, B)
%-   CHEBOP minus.
%   A - B subtracts the operators of the CHEBOPs A and B
%
%   C = MINUS(A, B) is called for the syntax 'A - B'.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

C = plus(A, uminus(B));
    
end
