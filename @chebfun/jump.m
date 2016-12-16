function j = jump(f, x, c)
%JUMP   The jump in a CHEBFUN over a breakpoint.
%   J = JUMP(F, X, C) is simply a wrapper for F(X, 'right') - F(X, 'left') - C.
%   If only two inputs are given, C is assumed to be zero.
%
% Example:
%   x = chebfun(@(x) x);
%   j = jump(sign(x), 0) % returns j = 2
% 
% See also FEVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    c = 0; 
end

% JUMP() is just a wrapper for FEVAL() and MINUS():
j = feval(f, x, 'right') - feval(f, x, 'left') - c;

end
