function [values, maskNaN, maskInf] = extrapolate(values)
%EXTRAPOLATE  Extrapolate data values from values at Chebyshev points of 
%   1st kind.
%
%   EXTRAPOLATE(VALUES) uses barycentric interpolants to extrapolate out 
%   NaNs or Infs from the numeric data in VALUES.
%
%   [VALUESOUT, MASKNAN, MASKINF] = EXTRAPOLATE(VALUESIN) returns logical
%   vectors indicating when a NaN or Inf was encountered in rows of VALUESIN.
%
%   Note that if any column of a multivalued function returns to NaN or Inf,
%   then _all_ columns are extrapolated at the point. Thus MASKNAN and MASKINF
%   are always column vectors, even if VALUES is a matrix.
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

    

