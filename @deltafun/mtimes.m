function f = mtimes(f, c)
%*   Multiplication of DELTAFUN objects.
%   F*C or C*F multiplies a DELTAFUN F by a scalar or matrix C.
% 
% See also TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % DELTAFUN * [] = []
    
    f = []; 
    return
    
elseif ( ~isa(f, 'deltafun') )       % First input is not a DELTAFUN
    
    % C must be a DELTAFUN and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);

elseif ( isa(c, 'double') )         % DELTAFUN * double 
    
    % DOUBLE*DELTAFUN requires that the double to be scalar.
    if ( ~isscalar(c) )
        error('CHEBFUN:DELTAFUN:mtimes:size', ...
              'Inner matrix dimensions must agree.');
    end
    
    % Multiply c with the smooth part:
    f.funPart = f.funPart * c;    
    f.deltaMag = f.deltaMag * c;
    
elseif ( isa(c, 'deltafun') )        % DELTAFUN * DELTAFUN  
    
    error('CHEBFUN:DELTAFUN:mtimes:deltafunMtimesdeltafun', ...
          'Use .* to multiply DELTAFUN objects.');
    
else                                % DELTAFUN * ???
    
    error('CHEBFUN:DELTAFUN:mtimes:deltafunMtimesUnknown', ...
          'mtimes does not know how to multiply a DELTAFUN and a %s.', class(c));
end

end
