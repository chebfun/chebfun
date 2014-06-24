function C = diffmat(N, varargin)
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
k = 1;
dom = [-1 1];
disc = colloc2();
for j = 1:numel(varargin)
    v = varargin{j};
    if ( isnumeric(v) )
        if ( isscalar(v) )
            k = v;
        else
            dom = v;
        end
    elseif ( isa(v, 'function_handle') || ischar(v) || ...
        isa(v, 'chebDiscretization') )
        disc = v;
    else
        error('CHEBFUN:TRUNK:diffmat:unknown', ...
            'Unknown input of type %s.', class(v));
    end
end
% Ensure DISC is a sicretization:
if ( ischar(disc) )
    disc = str2func(disc);
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
    