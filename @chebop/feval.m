function out = feval(N, varargin)
%FEVAL   Evaluate the operator of the CHEBOP at a CHEBFUN or CHEBMATRIX.
%   OUT = FEVAL(N, U) for a CHEBFUN or CHEBMATRIX U applies the CHEBOP N to U,
%   i.e., it returns N(U). Here, N.OP should be of the form @(u) diff(u,2) + ...
%   If N.op is of the form @(x, u) diff(u,2) + ... then an x variable is
%   instantiated internally and included automatically, however this should not
%   be relied upion.
%
%   OUT = FEVAL(N, X, U) for the CHEBFUN X and CHEBFUN or CHEBMATRIX U applies
%   the CHEBOP N to X and  U, i.e., it returns N(X, U) where N.OP has the form
%   @(x, u) diff(u,2) + .... Here, X shouyld be the dependent variable on
%   N.DOMAIN.
%
%   OUT = FEVAL(N, X, U1, U2, ..., UM) for a CHEBFUN X and CHEBFUN or CHEMBATRIX
%   objects U1, ..., UM applies the CHEBOP N to the functions Uk; i.e., it
%   returns N(X, U1, U2, ..., UM) where N.OP has the form @(x, u1, u2, ..., um).
%   Note that for systems of equations, X _must_ be included in N.OP.
%
%   OUT = FEVAL(N, X, U) where U is a CHEBMATRIX of M entries and N.OP has the
%   form @(X, U1, U2, ..., UM) is equivalent FEVAL(N, X, U{1}, ..., U{M}).
%   Again, OUT = FEVAL(N, U) will also work in this situation, but should not be
%   relied upon.
%
% See also CHEBOP/SUBSREF, LINOP/MTIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We must expand CHEBMATRIX entries out to a cell for {:} to work below.
isChebMat = cellfun(@(u) isa(u, 'chebmatrix'), varargin);
if ( any(isChebMat) )
    args = {};
    for k = 1:numel(varargin)
        % Append variables from the kth input:
        if ( isChebMat(k) )
            args = [args , varargin{k}.blocks.']; %#ok<AGROW>
        else
            args = [args , varargin(k)];          %#ok<AGROW>
        end
    end
else
    args = varargin;
end

numberOfInputs = nargin(N);
if ( numel(args) == numberOfInputs - 1 )
    % Check if we need to include an x (independent variable):
    x = chebfun(@(x) x, N.domain);
    args = [ {x} , args ];
elseif ( numberOfInputs == 0 )
    % Return empty:
    out = [];
    return
end

% Evaluate the operator:
out = N.op(args{:});

end
