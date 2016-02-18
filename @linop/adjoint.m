function L = adjoint(L, bctype)
%ADJOINT   Compute the adjoint of a LINOP.
%   L = ADJOINT(L), where L is a LINOP, returns the adjoint LINOP of L under
%   the assumption that L only has endpoint functional constraints.
%
%   L = ADJOINT(L, BCTYPE) allows for more general boundary conditions.
%
% See also ?.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check BCTYPE
if ( nargin < 2 )
    bctype = 'endpoints';
    warning('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    ['ADJOINT only supports purely differential operators with ',...
     'end point boundary conditions.']);


elseif ( ~strcmp(bctype,'periodic') && ~strcmp(bctype,'endpoints') )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'ADJOINT doesn''t support this type of boundary conditions for the moment.');

end

% Check diffOrder
if ( L.diffOrder < 0 )
    error('CHEBFUN:LINOP:adjoint:difforder', ...
    'ADJOINT doesn''t support integral operators for the moment.');

end

% Trivial case:
if ( L.diffOrder == 0 )
    return
end

% Formal adjoint
[Lstar, adjcoeffs] = formalAdjoint(L,bctype);

% Create adjoint linop
% periodic bcs
if ( strcmp(bctype,'periodic') )
    Lstar.constraint = L.constraint;
    L = Lstar;

% other bcs
else
    % Adjoint boundary conditions

end

end







function [L, adjcoeffs] = formalAdjoint(L,bctype)
%FORMALADJOINT   Computes the formal adjoint of a LINOP.
%   L = FORMALADJOINT(L, BCTYPE) returns the formal adjoint of a LINOP L.
%   If BCTYPE == 'periodic' then the adjoint uses periodic coefficients.
%
% See also adjointBCs.

% Get the value of the highest derivative:
n = L.diffOrder;

% Extract the blocks of the LINOP:
blocks = L.blocks;

% ADJOINT doesn't support system of equations for the moment.
% [TODO]: Support systems of equations.
if ( max(size(blocks,2)) > 1 )
    error('CHEBFUN:LINOP:adjoint:system', ...
        'ADJOINT doesn''t support systems of equations for the moment.');
end

% Get the first block and the domain:
block = blocks{1};
dom = L.domain;

% ADJOINT doesn't support piecewise domains at the moment.
% [TODO]: Support piecewise domains.
if ( size(dom,2) > 2 )
    error('CHEBFUN:LINOP:adjoint:domain', ...
        'ADJOINT doesn''t support piecewise domains for the moment.');
end

% Get the coefficients:
if ( n == 0 )
    % This is a multiplication operator, nothing to do here:
    if ( nargout > 1 )
        adjcoeffs = conj(block);
    end
    return
else
    coeffs = toCoeff(block);
end

% Compute the coefficients of the adjoint:
adjcoeffs = cell(n+1,1);
for k = 0:n
    adjcoeffs{k+1} = 0*coeffs{k+1};
end
for k = 0:n
    for l = 0:k
        adjcoeffs{n+1-l} = adjcoeffs{n+1-l} + ...
            (-1)^k*nchoosek(k,l)*conj(diff(coeffs{n+1-k}, k-l));
    end
end

% if periodic 
if ( strcmp(bctype,'periodic') )
    for k = 0:n
        c = adjcoeffs{k+1};
        adjcoeffs{k+1} = chebfun( @(x) feval(c,x), domain(c), 'trig');
    end
end

% Construct a LINOP from these new coefficients:
L = operatorBlock.mult(adjcoeffs{n},dom);
n = length(coeffs);
for k = 1:n-1
    L = L + operatorBlock.mult(adjcoeffs{n-k},dom)*operatorBlock.diff(dom,k);
end
L = linop(L);

if ( nargout > 1 )
    % Store the coefficients in a CHEBMATRIX:
    adjcoeffs = chebmatrix(adjcoeffs);
end

end

