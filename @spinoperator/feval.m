function out = feval(S, varargin)
%FEVAL   Evaluate the operator of a SPINOP (1, 2 or 3) at a CHEBFUN (1, 2, 3) or 
%at a CHEBMATRIX.
%   OUT = FEVAL(S, U) for a CHEBFUN (1, 2, or 3) or CHEBMATRIX U applies S to U.
%
%   OUT = FEVAL(S, U, V, ...) for CHEBFUN'S (1, 2, 3) U, V, ... applies S to
%   the CHEBFUN'S.
%
% See also SPINOPERATOR/SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

L = S.lin;
N = S.nonlin;

% CASE 1: S(U).
if ( nargin == 2 )
    
    u = varargin{1};
    
    if ( isa(u, 'chebfun') || isa(u, 'chebfun2') || isa(u, 'chebfun3') )
        
        out = L(u) + N(u);
        
    elseif ( isa(u, 'chebmatrix' ) )
        
        u = u.blocks.';
        out = L(u{:}) + N(u{:});
        
    end

% CASE 2: S(U,V,...).
else

    out = L(varargin{:}) + N(varargin{:});
    
end

% Convert to a CHEBMATRIX if it's a CHEBFUN2V or a CHEBFUN3V:
if ( isa(out, 'chebfun2v') || isa(out, 'chebfun3v') )
    
    u = chebmatrix(out(1));
    for k = 2:size(out, 1)
        u(k,1) = out(k);
    end
    out = u;
    
end

end