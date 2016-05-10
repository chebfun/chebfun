function varargout = deal(X, simplifyFlag)
%DEAL   Deal for chebmatrices.
%   [A,B,C,...] = DEAL(X) assigns to the variables A,B,C... the entries in the
%   chebmatrix X by row-major ordering so that A=X(1,1), B=X(1,2), C=X(1,3),
%   etc. If the number of outputs is less than the number of entries in X,
%   only those entries are assigned. By default, outputs are simplified (i.e.
%   they are passed through chebfun.simplify) before they are returned.
%
%   [A,B,C,...] = DEAL(X, SIMPLIFYFLAG), with SIMPLIFYFLAG = FALSE skips the
%   default simplification step.
%
%   DEAL is especially useful for the solution of differential equations,
%   wherein one may write, for example,
%
%       N = chebop(@(x,u,v) [diff(u) + v; diff(v) - u], [0,1]);
%       N.bc = @(x,u,v) [u(0) - 1; v(0)];
%       [u,v] = deal(N \ 0);  % solution components of a system
%                             % of differential equations
%
%   However, users do not normally need this syntax as it is
%   simpler to write  [u,v] = N \ 0 ; .


% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    % By default, call simplify:
    simplifyFlag = true;
end

if ( nargout > prod(size(X)) )
    error('CHEBFUN:CHEBMATRIX:deal:tooManyOutputs', ...
        ['Number of outputs should be less than or equal ' ...
         'to the number of elements in input chebmatrix.'])
end

X = X.';
blocks = X.blocks(:);

if ( simplifyFlag )
    for blockCounter = 1:length(blocks)
        if ( isa(blocks{blockCounter}, 'chebfun') )
            blocks{blockCounter} = simplify(blocks{blockCounter});
        end
    end
end

varargout = blocks(1:nargout);
