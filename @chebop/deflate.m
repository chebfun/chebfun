function N = deflate(N, r, p, alp)
%DEFLATE   Deflate known solutions from the operator of a CHEBOP
%
%   Calling sequence:
%       N = DEFLATE(N, R, P, ALP)
%   where the inputs are
%       N:    A CHEBOP, whose arguments are X and U
%       R:    A CHEBFUN or CHEBMATRIX of previously found solutions
%       P:    The power coefficient of the deflation scheme
%       ALP:  The shift coefficient of the deflation scheme
%   and the output is
%       N:    A CHEBOP, whose N.OP has had the solutions R deflated from
%
%   Note: DEFLATE(N, ...) only returns a modified CHEBOP, it is then necessary
%   to \ to compute new solutions.
%
%   Example (Painleve I):
%       L = 10;
%       d = [0 L];
%       N = chebop(@(x,u) diff(u,2)-u.^2+x, d);
%       N.lbc = 0;
%       N.rbc = sqrt(L);
%       r0 = N\0; % First solution
%       Ndef = deflate(N, r0, 3, .1);   % Deflate
%       r1 = Ndef\0;                    % Compute second solutions

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if (isa(r, 'chebfun'))
    % Ensure R is a CHEBMATRIX to pass to DEFLATIONFUN below:
    r = chebmatrix(r);
end

if (nargin(N.op) < 2 )
    error('must pass x and u to N.op');
end

N.op = @(x,u) deflationFun(N.op(x,u), u, r, p, alp);

end