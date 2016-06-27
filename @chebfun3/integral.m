function I = integral(f, varargin)
%INTEGRAL   Line integral of a CHEBFUN3 over a parametric curve.
%   I = INTEGRAL(F, G), returns the integral of the CHEBFUN3 object F
%   along a parametric curve defined by the Inf x 3 quasimatrix G. Columns 
%   of G represent parametrization of the 3D curve.
%
%   I = INTEGRAL(F), returns the triple definite integral of the CHEBFUN3 
%   object F over its domain of definition.
% 
% See also CHEBFUN3/INTEGRAL2, CHEBFUN3/INTEGRAL3, CHEBFUN3/SUM, 
% CHEBFUN3/SUM2, and CHEBFUN3/SUM3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) ) 
    I = [];
    return
end

if ( nargin == 1 )                       % Another way to do sum3(f) 
    % Triple definite integral:
    I = sum3(f); 
    
else
    if ( (isa(varargin{1}, 'chebfun')) && size(varargin{1}, 2) == 3 )  
        % Line integral over an Inf x 3 quasimatrix.
        % Get curve: 
        curve = varargin{1}; 
        xCurve = curve(:, 1); 
        yCurve = curve(:, 2); 
        zCurve = curve(:, 3);         
        diffC = diff(curve);
        ds_squared = diffC(:, 1).^2 + diffC(:, 2).^2 + diffC(:, 3).^2;
        I = sum(feval(f, xCurve, yCurve, zCurve) .* sqrt(ds_squared), ...
            curve.domain);
    else
        error('CHEBFUN:CHEBFUN3:integral:badInputs',...
        'Unrecognised input arguments.');
    end
end

end