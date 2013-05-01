function f = mtimes(f, c)
%*   Multiplication of BNDFUN objects.
%   F*C or C*F multiplies a BNDFUN F by a scalar or matrix C.
%
%   If F is an array-valued BNDFUN and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % BNDFUN * [] = []
    f = []; 
    return
elseif ( ~isa(f, 'bndfun') )      % BNDFUN is not the first input
    % DOUBLE*BNDFUN requires that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:BNDFUN:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a BNDFUN and F a scalar double. Call MTIMES again.
    c.onefun = mtimes(c.onefun, f);
    
    f = c;
    
elseif ( isa(c, 'double') )         % BNDFUN * double
    f.onefun = mtimes(f.onefun, c);
    
elseif ( isa(c, 'bndfun') )       % BNDFUN * BNDFUN
    error('CHEBFUN:BNDFUN:mtimes:BndfunMtimesBndfun', ...
        'Use .* to multiply BNDFUN objects.');
else                                % BNDFUN * ???
    error('CHEBFUN:BNDFUN:mtimes:bndfunMtimesUnknown',...
        'mtimes does not know how to multiply a BNDFUN and a %s.', class(c));
end

end
