function h = kron(f, g, varargin)
%KRON   Kronecker/outer product of two chebfuns.
%
% H = KRON(F,G) where F and G are array-valued chebfuns constructs a
% chebfun2.  If size(F) = [Inf, K] and size(G) = [K, Inf] then H is a rank K
% chebfun2 such that
%
%  H(x,y) = F(y,1)G(x,1) + ... + F(y,K)G(x,K).
%
% If size(F) = [K,Inf] and size(G) = [Inf, K] then H is a chebfun2 such that
%
%  H(x,y) = G(y,1)F(x,1) + ... + G(y,K)F(x,K).
%
% This is function analogue of the Matlab command KRON.
%
% See also KRON.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% TODO: Add a third argument to this command to allow the outerproduct of
% linops/chebops. 

if ( ~isa( f, 'chebfun' ) || ~isa( g, 'chebfun' ) )
    error('CHEBFUN:KRON:INPUTS','Both inputs should be chebfuns');
end

if ( isempty( f ) || isempty( g ) )
    h = chebfun2();  % return empty chebfun2
    return
end

% get domains of chebfuns
fint = f.domain;
gint = g.domain;

if ( length(fint) > 2 || length(gint) > 2 )
    error('CHEBFUN:KRON:BREAKPTS',...
        'The two chebfuns must be smooth and contain no break points.');
end

% check if we have the right sizes:
[mf, nf] = size( f );
[mg, ng] = size( g );

if ( ( mf ~= ng ) || ( nf ~=mg ) )
    error('CHEBFUN:KRON:SIZES',...
        'Inconsistent sizes for the Kronecker product of chebfuns.');
end

if ( isinf( nf ) )  % size(F) = [K,Inf] and size(G) = [Inf, K]
    h = chebfun2.outerProduct(g, f);
else                % size(F) = [Inf, K] and size(G) = [K, Inf]
    h = chebfun2.outerProduct(f, g);
end


end

