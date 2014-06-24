function C = diffmat(N, k, dom, disc)
%DIFFMAT   Differentiation matrix.
%   D = DIFFMAT(N) returns the NxN differentiation matrix associated with the
%   Chebyshev spectral collocation method at second-kind Chebyshev points. D =
%   DIFFMAT(N, K) returns the differentiation matrix of order K. See
%   COLLOC2.DIFFMAT for further details.
%
%   D = DIFFMAT(N, K, DOM) scales the differentiation matrix D to the domain
%   DOM. DOM should be a 1x2 vector.
%
%   D = DIFF(N, K, DOM, DISC) or DIFF(N, K, DISC) returns the differentiation
%   matrix associated with the tech DISC.
%
% See also DIFF, COLLOC2.DIFFMAT, CUMSUMMAT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Parse the inputs:
if ( (nargin < 2) || isempty(k))
    k = 1;
end
if ( (nargin == 3) )
    if ( isa(dom, 'function_handle') || ischar(dom) )
        disc = dom;
        dom = cheboppref().domain;
    else
        disc = colloc2();
    end
elseif ( nargin <= 2 )
    disc = colloc2();
    dom = cheboppref().domain;
end
% Ensure DISC is a sicretization:
if ( ischar(disc) )
    disc = str2fun(disc);
end
if ( isa(disc, 'function_handle') )
    disc = disc();
end
% No breakpoints allowed:
if ( numel(dom) > 2 )
    dom = dom([1 end]);
    warning('CHEBFUN:TRUNK:diffmat:noBreaks', ...
        'DIFFMAT does not support domains with  breakpoints.');
end

%% Call DISC.DIFFMAT(N) and scale appropriately.
scl = (2/(dom(end)-dom(1)))^k;
C = scl*disc.diffmat(N, k);

end
    