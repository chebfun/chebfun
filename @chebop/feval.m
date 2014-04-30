function Narg = feval(N, varargin)
%FEVAL Evaluate the operator of the CHEBOP at a CHEBFUN.
%
% The calling sequence to FEVAL for a CHEBOP N is determined by the calling
% sequence of N.OP. CHEBOP/FEVAL expects the same inputs as N.OP, in other
% words, if the independent variable X appears in the list of arguments for
% N.OP, CHEBOP/FEVAL requires that argument to be passed. Likewise, if X does
% not appear in the calling sequence to N.OP, X should not be passed in as an
% argument to CHEBOP/FEVAL.
%
% The method CHEBOP/FEVAL works as follows:
%
% OUT = FEVAL(N, U) for a CHEBFUN U or a CHEBMATRIX U applies the CHEBOP N to U,
% i.e., it returns N(U). Here, N.OP has to be of the form
%       N.op = @(u) diff(u,2) + ...
% or of the CHEBMATRIX form
%       N.op = @(u) [diff(u{1}) + u{2}; u{1} + diff(u{2})
%
% OUT = FEVAL(N, X, U) for the CHEBFUN X and a CHEBFUN U or a CHEBMATRIX U
% applies the CHEBOP N to X and  U, i.e., it returns N(X, U). Here, X is the
% dependent variable on the domain of A. Here, N.OP has to be of the form
%       N.op = @(x, u) diff(u,2) + ...
% or of the CHEBMATRIX form
%       N.op = @(x, u) [diff(u{1}) + u{2}; u{1} + diff(u{2})
% If N.OP is of this form, X has to be passed as an argument to FEVAL.
%
% OUT = FEVAL(N, X, U1, U2, ..., UM) for a CHEBFUN X and CHEBFUN objects
% U1, ..., UM applies the CHEBOP N to the functions Uk; i.e., it returns
% N(X, U1, U2, ..., UM). Here, X is the dependent variable on the domain of N.
% Here, N.OP has to be of the form
%       N.op = @(x, u, v) [diff(u) + v; u + v)
% If N.OP is of this form, X has to be passed as an argument to FEVAL (which is
% a general requirement of the CHEBOP class in any case).
%
% See also chebop/subsref, linop/mtimes.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Evaluate the operator!
Narg = N.op(varargin{:});

end
