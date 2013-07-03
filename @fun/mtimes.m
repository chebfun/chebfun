function f = mtimes(f, c)
%*   Multiplication of FUN objects.
%   F*C or C*F multiplies a FUN F by a scalar or matrix C.
%
%   If F is an array-valued FUN and C is a matrix of appropriate dimension, then
%   the natural matrix multiplication is performed.
%
%   Use .* (TIMES) to multiply two FUN objects.
%
% See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % FUN * [] = []
    
    % Return empty:
    f = []; 
    return
    
elseif ( ~isa(f, 'fun') )           % FUN is not the first input
    
    % DOUBLE*FUN requires that the double is scalar:
    if ( numel(f) > 1 )
        error('CHEBFUN:FUN:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a FUN and F a scalar double. Call ONEFUN/MTIMES().
    c.onefun = mtimes(c.onefun, f);
    
    % Assign c to be the output.
    f = c;
    
elseif ( isa(c, 'double') )         % FUN * double
    
    % Call MTIMES() at the ONEFUN level:
    f.onefun = mtimes(f.onefun, c);
    
elseif ( isa(c, 'fun') )            % FUN * FUN
    
    error('CHEBFUN:FUN:mtimes:FunMtimesFun', ...
        'Use .* to multiply FUN objects.');
    
else                                % FUN * ???
    
    error('CHEBFUN:FUN:mtimes:FunMtimesUnknown',...
        'mtimes does not know how to multiply a FUN and a %s.', class(c));
    
end

end