function h = kron(f, g)
%KRON   Kronecker/outer product of two chebfuns.
%   H = KRON(F,G) where F and G are array-valued CHEBFUN objects constructs a
%   CHEBFUN2.  If size(F) = [Inf, K] and size(G) = [K, Inf] then H is a rank K
%   CHEBFUN2 such that
%       H(x,y) = F(y,1)G(x,1) + ... + F(y,K)G(x,K).
%
%   If size(F) = [K,Inf] and size(G) = [Inf, K] then H is a chebfun2 such that
%       H(x,y) = G(y,1)F(x,1) + ... + G(y,K)F(x,K).
%
% See also KRON.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Add a third argument to this command to allow the outerproduct of
% LINOPs/CHEBOPs. 

if ( ( ~isa(f, 'chebfun') || ~isa(g, 'chebfun') ) )
    error('CHEBFUN:CHEBFUN:kron:inputs', 'Both inputs should be CHEBFUNs.');
end

if ( isempty(f) || isempty(g) )
    h = chebfun2();  % Return empty CHEBFUN2.
    return
end

% Get domains of CHEBFUNs:
fDom = f.domain;
gDom = g.domain;

if ( (length(fDom) > 2) || (length(gDom) > 2) )
    error('CHEBFUN:CHEBFUN:kron:breakpts',...
        'The two CHEBFUNs must be smooth and contain no break points.');
end

% check if we have the right sizes:
[mf, nf] = size( f );
[mg, ng] = size( g );

if ( (mf ~= ng) || (nf ~= mg) )
    error('CHEBFUN:CHEBFUN:kron:sizes',...
        'Inconsistent sizes for the Kronecker product of CHEBFUNs.');
end

if ( isinf( nf ) )  
    % size(F) = [K,Inf] and size(G) = [Inf, K]
    h = chebfun2.outerProduct(g, f);
    
else
    % size(F) = [Inf, K] and size(G) = [K, Inf]
    h = chebfun2.outerProduct(f, g);
    
end


end

