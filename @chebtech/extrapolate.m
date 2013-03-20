function [values, maskNaN, maskInf] = extrapolate(f)
%EXTRAPOLATE  Extrapolate data values from values at Chebyshev points of 
%   1st kind.
%
%   EXTRAPOLATE(F) uses barycentric interpolants to extrapolate out 
%   NaNs or Infs from the numeric data in F.VALUES.
%
%   [VALUES, MASKNAN, MASKINF] = EXTRAPOLATE(F) returns logical
%   vectors indicating when a NaN or Inf was encountered in rows of F.VALUES.
%
%   Note that if any column of a multivalued function returns NaN or Inf,
%   then _all_ columns are extrapolated at the point. Thus MASKNAN and MASKINF
%   are always column vectors, even if F.VALUES is a matrix.
%
%   The F.COEFFS field is not used/required; only F.VALUES is needed.
%
% See also PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

values = f.values;
maskNaN = any(isnan(values), 2);
maskInf = any(isinf(values), 2);
mask = maskNaN | maskInf;

if ( any(mask) )
    % Obtain Chebyshev points:
    n = size(values, 1);
    x = f.chebpts(n);

    % The good:
    xgood = x(~mask);
    if ( isempty(xgood) )
        error('CHEBFUN:CHEBTECH1:extrapolate:nans', ...
            'Too many NaNs to handle.')
    end
    
    % The bad:
    xnan = x(mask);
    
    % Compute the modified barycentric weights:
    w = f.barywts(n); % Standard weights.
    w = w(~mask);     % Barycentric weights corresponding to good points.
    for k = 1:length(xnan)
        % Compute the modified barycentric weights for the bad points:
        w = w.*( xgood - xnan(k) );
    end
    
    % Preallocate the storage for extrapolated values at the bad points:
    newvals = zeros(length(xnan), size(values, 2)); 
    % Barycentric formula of the second kind:
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
