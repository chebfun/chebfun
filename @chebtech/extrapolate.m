function [values, maskNaN, maskInf] = extrapolate(f, values)
%EXTRAPOLATE  Extrapolate data values from values at Chebyshev points of 
%             1st kind.
%   EXTRAPOLATE(F, VALUES) uses barycentric interpolants to extrapolate out 
%   NaNs or Infs from the numeric data in VALUES.
%
%   [VALUES, MASKNAN, MASKINF] = EXTRAPOLATE(F, VALUES) returns logical
%   vectors indicating when a NaN or Inf was encountered in rows of VALUES.
%
%   Note that if any column of an array-valued function returns NaN or Inf,
%   then _all_ columns are extrapolated at the point. Thus MASKNAN and MASKINF
%   are always column vectors, even if VALUES is a matrix.
%
%   The F.COEFFS field is not used/required; only VALUES is needed.
%
% See also PREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
        error('CHEBFUN:CHEBTECH:extrapolate:nansInfs', ...
            'Too many NaNs/Infs to handle.')
    end
    
    % The bad:
    xbad = x(mask);
    
    % Compute the modified barycentric weights:
    w = f.barywts(n); % Standard weights.
    w = w(~mask);     % Barycentric weights corresponding to the good points.
    for k = 1:length(xbad)
        % Compute the modified barycentric weights for the bad points:
        % (Recall that w_k = 1/prod_{j~=k)(x_j-x_k) )
        w = w.*( xgood - xbad(k) );
    end
    
    % Preallocate the storage for extrapolated values at the bad points:
    newvals = zeros(length(xbad), size(values, 2)); 
    % Barycentric formula of the second kind:
    for k = 1:length(xbad)
        % Compute the weights:
        w2 = w./(xbad(k) - xgood);
        % Sum the values:
        newvals(k,:) = (w2.'*values(~mask,:)) / sum(w2);
    end
    
    % Update the values:
    values(mask,:) = newvals;
    
end

end
