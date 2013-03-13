function [values, maskNaN, maskInf] = extrapolate(values)
%EXTRAPOLATE  Extrapolate data values from values at Chebyshev points.
%   EXTRAPOLATE(VALUES) uses barycentric interpolants to extrapolate out NaNs or
%   Infs from the numeric data in VALUES.
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

maskNaN = any(isnan(values), 2);
maskInf = any(isinf(values), 2);
mask = maskNaN | maskInf;

if ( any(mask) )
    
    % Obtain Chebyshev points
    n = size(values, 1);
    x = funcheb2.chebpts(n);
    
    if ( any(mask(2:end-1)) )   % Interior NaNs
        
        % Store information at the ends (will be replaced below if poss.)
        maskends = mask([1, end]);
        indxends = [1 ; n]; 
        indxends = indxends(~maskends);
        vends = values([1, end],:);
        
        % Force extrapolation at end points for now:
        mask([1, end]) = true;
    
        % The good:
        xgood = x(~mask);
        if ( isempty(xgood) )
            error('CHEBFUN:FUNCHEB2:extrapolate:nans', ...
                'Too many NaNs to handle.')
        end
        
        % The bad:
        xnan = x(mask);
        
        % Compute barycentric weights:
        w = ones(size(xgood));
        for k = 1:length(xnan)
            w = w.*abs(xnan(k) - xgood);
        end
        w(2:2:end) = -w(2:2:end);
        
        % Mini barycentric formula:
        newvals = zeros(length(xnan), size(values, 2));
        for k = 1:length(xnan)
            w2 = w./(xnan(k) - xgood);
            newvals(k,:) = (w2.'*values(~mask,:)) / sum(w2);
        end
        
        % Update the values:
        values(mask,:) = newvals;
        
        % Keep the ends if possible:
        values(indxends,:) = vends(~maskends,:);
        
    else               % Force extrapolation at endpoints anyway
        
        % Extrapolate at endpoints if needed using "Fejer's 2nd rule" type of
        % barycentric formula.
        
        xi = x(2:end-1); % Interior nodes
        
        % Barycentric weights for left evaluation.
        w = 1 - xi;                 
        w(2:2:end) = -w(2:2:end);
        values(1,:) = (w.'*values(2:end-1,:)) / sum(w);   % Values at x = -1;
        
        % Barycentric weights for right evaluation.
        w = 1 + xi;                 
        w(2:2:end) = -w(2:2:end);
        values(end,:) = (w.'*values(2:end-1,:)) / sum(w); % Values at x = 1;

        
    end

end

end

    

