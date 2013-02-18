function f = mtimes(f, c)
%*	Multiplication of FUNCHEB2 objects.
%   F*C or C*F multiplies a FUNCHEB2 F by a scalar or matrix C.
%
%   If F is a verctorised FUNCHEB2 and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % FUNCHEB2 * [] = []
    f = []; 
    return
elseif ( ~isa(f, 'funcheb2') )      % FUNCHEB2 is not the first input

    % DOUBLE*FUNCHEB2 will require that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:FUNCHEB2:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a FUNCHEB2 and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % FUNCHEB2 * double
    
    % Check dimensions
    if ( size(f.values, 2) ~= size(c, 1) && numel(c) > 1)
        error('CHEBFUN:FUNCHEB2:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    f.vscale = norm(c(:), inf)*f.vscale;
    
    % If the vertical scale is zero, then set the FUNCHEB2 to zero:
    if ( f.vscale == 0 )
        f.values = zeros(size(f.values, 1), 1);
        f.coeffs = zeros(size(f.values, 1), 1);
    end
    
elseif ( isa(c, 'funcheb2') )       % FUNCHEB2 * FUNCHEB2
    
    error('CHEBFUN:FUNCHEB2:mtimes:funcheb2MtimesFun2', ...
        'Use .* to multiply FUNCHEB2s.');
    
else                                % FUNCHEB2 * ???
    
    error('CHEBFUN:FUNCHEB2:mtimes:funcheb2MtimesUnknown',...
        'mtimes does not know how to multiply a FUNCHEB2 and a %s.', class(c));
    
end
