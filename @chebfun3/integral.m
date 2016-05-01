function I = integral(f, varargin)
%INTEGRAL   Line integral of a CHEBFUN3 over a curve.
%
%   I = INTEGRAL(f), returns the definite integral of a CHEBFUN3 over its 
%   domain of definition.
% 
%   I = INTEGRAL(f, G), returns the integral of the CHEBFUN3 object f
%   along the curve defined by the inf x 3 quasimatrix G. Columns of G 
%   represent parametrization of the 3D curve. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )                         % Another way to do sum3(f) 
    %I = integral3(f); 
    I = sum3(f); 
    
else
    if ( (isa(varargin{1}, 'chebfun')) && size(varargin{1}, 2) == 3 )  
        % Line integral over an inf x 3 quasimatrix
        % Get curve: 
        curve = varargin{1}; 
        xcurve = curve(:, 1); 
        ycurve = curve(:, 2); 
        zcurve = curve(:, 3);         
        diffC = diff(curve);
        ds_squared = diffC(:, 1).^2 + diffC(:, 2).^2 + diffC(:, 3).^2;
        I = sum(feval(f, xcurve, ycurve, zcurve) .* sqrt(ds_squared), ...
            curve.domain);
    else
        error('CHEBFUN:CHEBFUN3:integral:badInputs',...
        'Unrecognised input arguments.');
    end
end

end