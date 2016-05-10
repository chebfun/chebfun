function f = mtimes(f, c)
%*   Multiplication of SINGFUN objects.
%   F*C or C*F multiplies a SINGFUN F by a scalar C.
%
% See also TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % SINGFUN * [] = []
    f = []; 
    
elseif ( ~isa(f, 'singfun') )       % First input is not a SINGFUN
   
    % C must be a SINGFUN and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    
elseif ( isa(c, 'double') )         % SINGFUN * double  

    % SINGFUN*DOUBLE requires that the double is scalar.
    if ( numel(c) > 1 )
        error('CHEBFUN:SINGFUN:mtimes:size', ...
              'Inner matrix dimensions must agree.');
    end
    
    % Multiply c with the smooth part:
    f.smoothPart = f.smoothPart * c;    
    
elseif ( isa(c, 'singfun') )        % SINGFUN * SINGFUN  
    error('CHEBFUN:SINGFUN:mtimes:singfunMtimesSingfun', ...
          'Use .* to multiply SINGFUN objects.');
    
else                                % SINGFUN * ???
    error('CHEBFUN:SINGFUN:mtimes:singfunMtimesUnknown', ...
          'mtimes does not know how to multiply a SINGFUN and a %s.', class(c));
end

end
