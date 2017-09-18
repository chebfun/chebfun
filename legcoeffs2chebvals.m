function v_cheb = legcoeffs2chebvals(c_leg, varargin)
%LEGCOEFFS2CHEBVALS  Convert Legendre coefficients to Chebyshev values. 
%   V_CHEB = LEGCOEFFS2CHEBVALS(C_LEG) converts the vector C_LEG of Legendre
%   coefficients to a vector V_CHEB of values at second-kind Chebyshev points.
%
%   V_CHEB = LEGCOEFFS2CHEBVALS(C_LEG, 1) is similar, but evaluates at
%   first-kind Chebyshev points.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
v1 = {}; v2 = {};
for k = 1:numel(varargin)
    if ( isnumeric(varargin{k}) )
        v1{1} = varargin{k};
    elseif ( ischar(varargin{k}) )
        v2{1} = varargin{k};
    end
end

c_cheb = leg2cheb(c_leg, v2{:}); % Convert to Chebyshev coefficients.
v_cheb = chebcoeffs2chebvals(c_cheb, v1{:});  % Convert to Chebyshev values.

end
