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
% See also funcheb.pref.

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
    
    % Compute the modified barycentric weights:
    w = funcheb1.barywts(n); % Standard weights.
    w = w(~mask); % Barycentric weights corresponding to good points.
    for k = 1:length(xnan)
        % Compute the modified barycentric weights for the bad points:
        w = w.*( xgood - xnan(k) );
    end
    
    % Preallocate the storage for extrapolated values at the bad points.
    newvals = zeros(length(xnan), size(values, 2)); 
    % Barycentric formula of the second (true) kind:
    for k = 1:length(xnan)
        % Compute the weights:
        w2 = w./(xnan(k) - xgood);
        % Sum the values:
        newvals(k,:) = (w2.'*values(~mask,:)) / sum(w2);
    end
    
    % Update the values:
    values(mask,:) = newvals;
    
end

end

    

