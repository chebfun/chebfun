function f = mtimes(f, c)
%*   Multiplication of CLASSICFUN objects.
%   F*C or C*F multiplies a CLASSICFUN F by a scalar or matrix C.
%
%   If F is an array-valued CLASSICFUN and C is a matrix of appropriate dimension, then
%   the natural matrix multiplication is performed.
%
%   Use .* (TIMES) to multiply two CLASSICFUN objects.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % CLASSICFUN * [] = []
    
    % Return empty:
    f = []; 
    return
    
elseif ( ~isa(f, 'classicfun') )           % CLASSICFUN is not the first input
    
    % DOUBLE*CLASSICFUN requires that the double is scalar:
    if ( numel(f) > 1 )
        error('CHEBFUN:CLASSICFUN:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a CLASSICFUN and F a scalar double. Call ONEFUN/MTIMES().
    c.onefun = mtimes(c.onefun, f);
    
    % Assign c to be the output.
    f = c;
    
elseif ( isa(c, 'double') )         % CLASSICFUN * double
    
    % Call MTIMES() at the ONEFUN level:
    f.onefun = mtimes(f.onefun, c);
    
elseif ( isa(c, 'classicfun') )            % CLASSICFUN * CLASSICFUN
    
    error('CHEBFUN:CLASSICFUN:mtimes:classicfunMtimesClassicfun', ...
        'Use .* to multiply CLASSICFUN objects.');
    
else                                % CLASSICFUN * ???
    
    error('CHEBFUN:CLASSICFUN:mtimes:classicfunMtimesUnknown',...
        'mtimes does not know how to multiply a CLASSICFUN and a %s.', class(c));
    
end

end
