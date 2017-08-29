function c_leg = chebvals2legcoeffs(v_cheb, varargin)
%CHEBVALS2LEGCOEFFS  Convert Chebyshev values to Legendre coefficients.
% 	LEGCOEFFS = CHEBVALS2LEGCOEFFS(CHEBVALS), converts the column vector
% 	CHEBVALS of values on a second-kind Chebyshev grid (i.e, F(CHEBPTS(N))) to a
% 	vector LEGCOEFFS of Legendre coefficients, where the degree k Legendre
% 	polynomial P{k} is normalized so that max(|P{k}|) = 1.
%
%   CHEBVALS2LEGCOEFFS(CHEBVALS, 1) is similar, but assumes the entries in
% 	CHEBVALS come from evaluating on a first-kind Chebyshev grid, i.e., 
%   F(CHEBPTS(N, 1))).
%
%   CHEBVALS2LEGCOEFFS(CHEBVALS, 'norm') or CHEBVALS2LEGCOEFFS(CHEBVALS, 1,
%   'norm')  is as above, but with the Legendre polynomials normalized to be
%   orthonormal.
%
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

c_cheb = chebvals2chebcoeffs(v_cheb, v1{:}); % Convert to Chebyshev coefficients.
c_leg = cheb2leg(c_cheb, v2{:});             % Convert to Legendre coefficients.

end
