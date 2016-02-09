function N = deflate(N, r, p, alp, type)
%DEFLATE   Deflate known solutions from the operator of a CHEBOP.
%
%   Calling sequence:
%       N = DEFLATE(N, R, P, ALP)
%   where the inputs are
%       N:    A CHEBOP, whose arguments are X and U. At the moment, only scalar
%             problems are supported.
%       R:    A CHEBFUN or CHEBMATRIX of previously found solutions
%       P:    The power coefficient of the deflation scheme
%       ALP:  The shift coefficient of the deflation scheme
%       TYPE: Type of norm used for deflation. Possible options are 'L2'
%             (default) and 'H1'.
%   and the output is
%       N:    A CHEBOP, whose N.OP has had the solutions R deflated from
%
%   The deflation operator constructed corresponds to a new CHEBOP G, so that
%
%       G(u) = M(u, r)N(u)
%   
%   where N is the input CHEBOP, and M is the deflation operator
%       
%       M(u;r) = 1/(||u-r||^p_{type}) + alp
%
%   or in the case of multiple roots r_1, ..., r_n being deflated
%
%       M(u; r_1, ... r_n) = 1/(||u-r_1||^p_{type}*...*||u-r_n||^p_{type}) + alp
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
%
% References:
%   [1] Deflation techniques for finding distinct solutions of nonlinear
%   partial differential equations (P. E. Farrell, Á. Birkisson, S. W. Funke),
%   In SIAM Journal on Scientific Computing, volume 37, 2015.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 5 )
    type = 'L2';
end

if (isa(r, 'chebfun'))
    % Ensure R is a CHEBMATRIX to pass to DEFLATIONFUN below:
    r = chebmatrix(r);
end

assert(nargin(N.op) <= 2, 'CHEBFUN:CHEBOP:deflate:nargin', ...
    'Currently, CHEBOP deflation only supports scalar problems.');

% Modify the output operator to reflect the deflation:
if (nargin(N.op) == 1 )
    N.op = @(x,u) deflationFun(N.op(u), u, r, p, alp, type);
else
    N.op = @(x,u) deflationFun(N.op(x,u), u, r, p, alp, type);
end

end
