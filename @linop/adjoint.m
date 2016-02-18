function [Lstar, adjcoeffs, bcOp] = adjoint(L, bcType)
%ADJOINT   Compute the adjoint of a LINOP.
%   ADJOINT(L), where L is a LINOP, returns the adjoint LINOP of L under
%   the assumption that L only has endpoint or periodic functional constraints.
%
%   ADJOINT(L, BCTYPE) allows for more general boundary conditions.
%
% See also ?.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( L.diffOrder == 0 )
    Lstar = L;
    return
end

%% 
% Check to so if we can do compute an adjoint for this operator!
parseInputs(L, bcType);

%%
% Formal adjoint:
pref = chebfunpref();
if ( strcmpi(bcType, 'periodic') )
    pref.tech = @trigtech;
end 
[Lstar, adjcoeffs] = formalAdjoint(L, pref);

%%
% Create adjoint linop
if ( strcmp(bcType, 'periodic') )
    % Periodic bcs
    Lstar.constraint = L.constraint;
    B = 'periodic';
else
    % Adjoint boundary conditions
    [constraint, B] = adjointBCs(L, bcType);
    Lstar.constraint = constraint;
end

if ( nargout == 3 )
    bcOp = prettyPrint(B);
end

end

function parseInputs(L, bcType)

if ( L.diffOrder < 0 )
    % Check for integral operators:
    error('CHEBFUN:LINOP:adjoint:difforder', ...
    'ADJOINT doesn''t support integral operators for the moment.')
elseif ( size(L.domain, 2) > 2 )
    % [TODO]: Support piecewise domains.
    error('CHEBFUN:LINOP:adjoint:domain', ...
        'ADJOINT doesn''t support piecewise domains for the moment.');
elseif ( ~any(strcmp(bcType, {'periodic', 'bvp'})) )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'ADJOINT doesn''t support this type of boundary conditions for the moment.');
elseif ( all(size(L.blocks)) ~= 1 )
    % [TODO]: Support systems of equations.
    error('CHEBFUN:LINOP:adjoint:notscalar', ...
    'ADJOINT only support square scalar systems..');
end

end

function [L, adjCoeffs] = formalAdjoint(L, pref)
%FORMALADJOINT   Computes the formal adjoint of a LINOP.
%   L = FORMALADJOINT(L, BCTYPE) returns the formal adjoint of a LINOP L.
%
% See also adjointBCs.

% Get the domain and the value of the highest derivative:
dom = L.domain;
n = L.diffOrder;

% Get the coefficients:
coeffs = toCoeff(L.blocks{1}, pref);

% Compute the coefficients of the adjoint:
adjCoeffs = 0*coeffs;
for k = 0:n
    for l = 0:k
        adjCoeffs(n+1-l) = adjCoeffs{n+1-l} + ...
            (-1)^k*nchoosek(k,l)*conj(diff(coeffs{n+1-k}, k-l));
    end
end

% Construct a LINOP from these new coefficients:
L = 0;
M = @(f) operatorBlock.mult(f, dom);
D = @(k) operatorBlock.diff(dom, k);
for k = 0:n
    L = L + M(adjCoeffs{n+1-k}) * D(k);
end
L = linop(L);

end

function [constraint, B] = adjointBCs(L, bcType)

end

function [bcOpL, bcOpR] = prettyPrint(B)

bcOpL = '';
bcOpR = '';

end
