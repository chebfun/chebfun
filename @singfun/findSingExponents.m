function exponents = findSingExponents(op, singType)
%FINDSINGEXPONENTS   Endpoint singularity detection by sampling values.
%   Private method of SINGFUN.
%
%   This method is a wrapper for the underlying algorithms used to detect
%   orders of singularities in a SINGFUN.
%
% See also FINDPOLEORDER, FINDSINGORDER.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
exponents = zeros(1, 2);

if ( isnumeric(op) )
    return
end

% loop through each end
singEnd = {'left', 'right'};
for k = 1:2
    if ( strcmpi(singType{k}, 'pole') )
        exponents(k) = singfun.findPoleOrder(op, singEnd{k});
        
    elseif ( any(strcmpi(singType{k}, {'sing', 'root'})) )
        exponents(k) = singfun.findSingOrder(op, singEnd{k});
        
    elseif ( strcmpi(singType{k}, 'none') )
        exponents(k) = 0;
        
    else
        error('CHEBFUN:SINGFUN:findSingExponents:unknownPref', ...
              'singType "%s" unknown.', singType{k})
          
    end    
end

end
