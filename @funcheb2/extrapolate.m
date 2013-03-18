function [values, maskNaN, maskInf] = extrapolate(values)
%EXTRAPOLATE  Extrapolate data values from values at Chebyshev points.
%   EXTRAPOLATE(VALUES) uses barycentric interpolation to extrapolate out NaNs 
%   or Infs from the numeric data in VALUES.
%
%   [VALUESOUT, MASKNAN, MASKINF] = EXTRAPOLATE(VALUESIN) returns logical
%   vectors indicating when a NaN or Inf was encountered in rows of VALUESIN.
%
%   Note that if any column of a multivalued function returns to NaN or Inf,
%   then _all_ columns are extrapolated at the point. Thus MASKNAN and MASKINF
%   are always column vectors, even if VALUES is a matrix.
%
% See also FUNCHEB.pref.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebun.org/ for Chebfun information.

% Bookkeep all bad points including NaNs and Infs.

maskNaN = any(isnan(values), 2);
maskInf = any(isinf(values), 2);
mask = maskNaN | maskInf;

if ( any(mask) ) % Do extrapolation if there is any bad point.
    
    % Compute the Chebyshev points of 2nd kind
    n = size(values, 1);
    x = funcheb2.chebpts(n);
    
    % If there is any interior bad point, do extrapolation at BOTH interior
    % point(s) and endpoints (i.e. -1 and 1) using the barycentric weights 
    % computed by calling barywts@funcheb2. Otherwise only the value of the
    % function at the endpoints (i.e. -1 and 1) is extrapolated using newly
    % computed barycentric weights, since the weights can be much more 
    % easily and simply obtained in this way, instead of by calling barywts.    
    
    if ( any(mask(2:end-1)) )   % Interior bad points, if there is any.
        
        % Force extrapolation at the endpoints for time-being.
        mask([1, end]) = true;
    
        % The good points:
        xgood = x(~mask);
        if ( isempty(xgood) )
            error('CHEBFUN:FUNCHEB2:extrapolate:nans', ...
                'Too many NaNs to handle.')
        end
        
        % The bad points:
        xnan = x(mask);
        
        % Compute the modified barycentric weights:
        w = funcheb2.barywts(n); % barycentric weights
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
        
    else % Force extrapolation at the endpoints anyway, even though there isn't any interior bad point. 
        
        % We compute the modified barycentric weights for extrapolating the
        % endpoints from the scratch instead of calling barywts.
        
        xi = x(2:end-1); % Interior nodes
        
        % Barycentric weights for extrapolating at -1.
        w = 1 - xi;                 
        w(2:2:end) = -w(1:2:end);
        values(1,:) = (w.'*values(2:end-1,:)) / sum(w);   % Values at x = -1;
        
        % Barycentric weights for extrapolating at 1.
        w = 1 + xi;                 
        w(2:2:end) = -w(1:2:end);
        values(end,:) = (w.'*values(2:end-1,:)) / sum(w); % Values at x = 1;
        
    end

end

end

    

