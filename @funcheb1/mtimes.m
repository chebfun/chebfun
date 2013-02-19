function f = mtimes(f, c)
%*	Scalar multiplication of FUNCHEB1 objects.
%   F*C or C*F multiplies a FUNCHEB1 F by a scalar C.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% FUNCHEB1 * [] = []
if ( isempty(f) || isempty(c) )
    f = []; 
    return
end

if ( ~isa(f, 'funcheb1') )      % FUNCHEB1 is not the first input

    % DOUBLE*FUNCHEB1 will require that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:FUNCHEB1:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a FUNCHEB1 and F a scalar double. Call MTIMES again.
    f = mtimes(c,f);
    return
    
elseif ( isa(c, 'double') ) % FUNCHEB1 * double
    
    % Check dimensions
    if ( size(f.values, 2) ~= size(c, 1) )
        error('CHEBFUN:FUNCHEB1:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    f.vscale = norm(c, inf)*f.vscale;
    
    % If the vertical scale is zero, then set the FUNCHEB1 to zero:
    if ( f.vscale == 0 )
        f.values = 0;
        f.coeffs = 0;
    end
    
elseif ( isa(c, 'funcheb1') )    % FUNCHEB1 * FUNCHEB1
    
    error('CHEBFUN:FUNCHEB1:mtimes:funcheb1Mtimesfuncheb1', 'Use .* to multiply FUNCHEB1s.');
    
else                        % FUNCHEB1 * ???
    
    error('CHEBFUN:FUNCHEB1:mtimes:funcheb1MtimesUnknown',...
        'mtimes does not know how to multiply a FUNCHEB1 and a %s.', class(c));
    
end
