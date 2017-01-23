function I = integral(f, varargin)
%INTEGRAL   Complete definite integral of SEPARABLEAPPROX.
%
%   I = INTEGRAL(F), returns the definite integral of a SEPARABLEAPPROX
%   integrated over its domain of definition.
%
%   I = INTEGRAL(F, g), returns the integral of a SEPARABLEAPPROX along the
%   curve defined by the real- or complex-valued CHEBFUN g.
%
% See also INTEGRAL2, SUM2, QUAD2D.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )                         % Another way to do sum2(f)
    
    I = integral2(f);
    
else
    
    if ( ~isa(varargin{1}, 'chebfun') )  % Line integral over a CHEBFUN
        I = integral2(f, varargin{ : });
    else
        % Get curve:
        c = varargin{1};
        % Line integral:
        I = sum(compose(c, f) .* abs(diff(c)));
    end
    
end

end
