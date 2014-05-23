function exponents = findSingExponents(op, singType, lr)
%FINDSINGEXPONENTS   Endpoint singularity detection by sampling values.
%   Private method of SINGFUN.
%
%   EXPONENTS = FINDSINGEXPONENTS(OP, SINGTYPE) tries to determine the singularity
%   strength at both endpoints.
%
%   EXPONENTS = FINDSINGEXPONENTS(OP, SINGTYPE, LR) determines the singularity
%   strength at the endpoint indicated by LR (1 for the left endpoint, 2 for the
%   right endpoint).
%
%   This method is a wrapper for the underlying algorithms used to detect
%   orders of singularities in a SINGFUN.
%
% See also FINDPOLEORDER, FINDSINGORDER.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If only two arguments are passed, try to determine singularity strength at
% both endpoints.
if ( nargin == 2 )
    singEnd = {'left', 'right'};
end

% If the third argument is given, only one endpoint needs to be determined.
if ( nargin > 2 )
    if ( lr == 1 )
        singEnd = {'left'};
    elseif ( lr == 2 )
        singEnd = {'right'};
    else
        error('CHEBFUN:SINGFUN:findSingExponents:wrongIdx', ...
              'The indicator for endpoints is not understandable.')
    end
end

for k = 1:numel(singEnd)
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
