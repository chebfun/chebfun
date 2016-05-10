function h = kron(f, g, varargin)
%KRON   Kronecker/outer product of two CHEBFUNs.
%
% H = KRON(F,G) where F and G are array-valued CHEBFUN objects constructs a
% CHEBFUN2.  If size(F) = [Inf, K] and size(G) = [K, Inf] then H is a rank K
% CHEBFUN2 such that
%
%   H(x,y) = F(y,1)G(x,1) + ... + F(y,K)G(x,K).
%
% If size(F) = [K,Inf] and size(G) = [Inf, K] then H is a chebfun2 such that
% 
%   H(x,y) = G(y,1)F(x,1) + ... + G(y,K)F(x,K).
%
% This is function analogue of the MATLAB command KRON.
%
% H = KRON(F, G, 'op') or H = KRON(F, G, 'operator') if size(F) = [Inf, K] and
% size(G) = [K, Inf], results in a rank-k OPERATORBLOCK A such that 
%   A*U = F*(G*U)
% for any CHEBFUN U.
%
% See also KRON.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

assert( isa(f, 'chebfun') && isa(g, 'chebfun'), ...
    'CHEBFUN:CHEBFUN:kron:inputs', ...
    'Both inputs should be CHEBFUNs.')

if ( isempty(f) || isempty(g) )
    h = chebfun2();  % Return empty CHEBFUN2.
    return
end

% Get domains of CHEBFUNs, and check that they don't contain breakpoints:
fDom = f.domain;
gDom = g.domain;
assert( length(fDom) <=2 && length(gDom) <=2, ...
    'CHEBFUN:CHEBFUN:kron:breakpts', ...
     'The two CHEBFUNs must be smooth and contain no break points.')

% Check if we have the right sizes:
[mf, nf] = size( f );
[mg, ng] = size( g );
assert( (mf == ng) && (nf == mg), ...
    'CHEBFUN:CHEBFUN:kron:sizes',...
    'Inconsistent sizes for the Kronecker product of CHEBFUNs.');

if ( nargin <= 2 )

    if ( isinf( nf ) )
        % size(F) = [K,Inf] and size(G) = [Inf, K]
        h = chebfun2.outerProduct(g, f);
        
    else
        % size(F) = [Inf, K] and size(G) = [K, Inf]
        h = chebfun2.outerProduct(f, g);
    end
    
elseif ( nargin == 3 )
    % Developers note: Historically, this used to be called by F*G', where F and
    % G where two CHEBFUNs.  Thus functionality has been moved to here since
    % F*G' now returns a low rank CHEBFUN2. Also, we haven't come across an
    % application for this method for quasimatrices (rather than array valued
    % chebfuns), which is why we allow us just to throw an error in that case.
    
    assert( strcmpi(varargin{1}, 'op') || strcmpi(varargin{1}, 'operator'), ...
        'CHEBFUN:CHEBFUN:kron:opts', ...
        'Unrecognized optional parameter.');

    assert( all(fDom == gDom), ...
        'CHEBFUN:CHEBFUN:kron:domain',...
        'Domains must be identical for Kronecker products returning operators.')

    assert( (numel(f) == 1) && (numel(g) == 1), ...
        'CHEBFUN:CHEBFUN:kron:quasimatrix', ...
        ['chebfun/kron with an operator output currently only supports ' ...
         'array-valued chebfuns, not quasimatrices.'])
     
    assert( isinf(mf) && isinf(ng) , ...
        'CHEBFUN:CHEBFUN:kron:columnsAndRows', ...
        ['chebfun/kron with an operator output only supports inputs such ' ...
         'that the first input is a column chebfun, and the second is a ' ...
         'row chebfun.'])
    
    % Construct the OPERATORBLOCK to be returned:
    h = operatorBlock.outer(f, g, fDom);

else
    error('CHEBFUN:CHEBFUN:kron:nargin','Too many input arguments.');
end
    

end

