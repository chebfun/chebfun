function C = cumsummat(N, dom, disc)
%CUMSUMMAT   Indefinite integration matrix.
%   C = CUMSUMMAT(N) returns the NxN indefinite integration matrix associated
%   with the Chebyshev spectral collocation method at second-kind Chebyshev
%   points. By convection, the arbitrary constant is chosen so that the result
%   is zero at -1. See CHEBCOLLOC2.CUMSUMMAT for further details.
%
%   D = CUMSUMMAT(N, DOM) scales the indefinite integration matrix D to the
%   domain DOM. DOM should be a 1x2 vector.
%
%   D = CUMSUMMAT(N, DOM, DISC) or CUMSUMMAT(N, DISC) returns the indefinite
%   integration matrix associated with the CHEBDISCRETIZATION DISC.
%
% See also CUMSUM, CHEBCOLLOC2.DIFFMAT, DIFFMAT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Parse the inputs:
if ( (nargin == 2) )
    if ( isa(dom, 'function_handle') || ischar(dom) )
        disc = dom;
        dom = cheboppref().domain;
    else
        disc = chebcolloc2();
    end
elseif ( nargin == 1 )
    disc = chebcolloc2();
    dom = cheboppref().domain;
end
% Ensure DISC is a discretization:
if ( ischar(disc) )
    disc = str2func(disc);
end
if ( isa(disc, 'function_handle') )
    disc = disc();
end
% No breakpoints allowed:
if ( numel(dom) > 2 )
    dom = dom([1 end]);
    warning('CHEBFUN:cumsummat:noBreaks', ...
        'CUMSUMMAT does not support domains with breakpoints.');
end

%% Call DISC.DIFFMAT(N) and scale appropriately.
scl = .5*(dom(end) - dom(1));
C = scl*disc.cumsummat(N);

end
    
