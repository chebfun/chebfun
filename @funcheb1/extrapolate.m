function [values, maskNaN, maskInf] = extrapolate(values)
%EXTRAPOLATE  Extrapolate data values from values at Chebyshev points.
%   EXTRAPOLATE(VALUES) extrapolates the given data VALUES at data points
%   wherever a NaN or an Inf is encountered using the barycentric interpolant
%   through the remaining numeric values. 
%
%   EXTRAPOLATE(VALUES, PREF) will do the same, as well as auotmatically
%   extrapolate the data at -1 and +1 if the pref.funcheb1('extrapolate', true).
%
%   [VALUESOUT, MASKNAN, MASKINF] = EXTRAPOLATE(VALUESIN) returns also logical
%   vectors indicating when NaNs or Infs were encountered in VALUESIN.
%
% See also FUNCHEB1.pref.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebun.org/ for Chebfun information.

maskNaN = any(isnan(values), 2);
maskInf = any(isinf(values), 2);
mask = maskNaN | maskInf;

if ( any(mask) )
    
    % Obtain Chebyshev points
    n = size(values, 1);
    x = funcheb1.chebpts(n);

    % The good:
    xgood = x(~mask);
    if ( isempty(xgood) )
        error('CHEBFUN:FUNCHEB1:extrapolate:nans', ...
            'Too many NaNs to handle.')
    end
    
    % The bad:
    xnan = x(mask);
    
    % Compute barycentric weights:
    w = funcheb1.barywts(n);
    w = w(~mask);
    for k = 1:length(xnan)
        w = w.*( xgood - xnan(k) );
    end
    
    % Barycentric formula of the second (true) kind:
    newvals = zeros(length(xnan), size(values, 2));
    for k = 1:length(xnan)
        w2 = w./(xnan(k) - xgood);
        newvals(k,:) = (w2.'*values(~mask,:)) / sum(w2);
    end
    
    % Update the values:
    values(mask,:) = newvals;
    
end

end

    

