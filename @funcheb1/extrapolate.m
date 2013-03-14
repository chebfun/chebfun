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

% Bookkeep all bad points including NaNs and Infs.

maskNaN = any(isnan(values), 2);
maskInf = any(isinf(values), 2);
mask = maskNaN | maskInf;

if ( any(mask) )  % Do extrapolation if there is any bad point.
    
    % Compute the Chebyshev points of 1st kind
    n = size(values, 1);
    x = funcheb1.chebpts(n);

    % The good points:
    xgood = x(~mask);
    if ( isempty(xgood) )
        error('CHEBFUN:FUNCHEB1:extrapolate:nans', ...
            'Too many NaNs to handle.')
    end
    
    % The bad points:
    xnan = x(mask);
    
    % Compute the modified barycentric weights:
    w = funcheb1.barywts(n); % barycentric weights
    w = w(~mask); % barycentric weights corresponding to good points
    for k = 1:length(xnan)
        w = w.*( xgood - xnan(k) ); % compute the modified barycentric weights
    end
    
    % Barycentric formula of the second (true) kind:
    newvals = zeros(length(xnan), size(values, 2)); % preallocate the storage for extrapolated values at the bad points.
    for k = 1:length(xnan)
        w2 = w./(xnan(k) - xgood); % compute the weights
        newvals(k,:) = (w2.'*values(~mask,:)) / sum(w2);
    end
    
    % Update the values at the bad points:
    values(mask,:) = newvals;
    
end

end

    

