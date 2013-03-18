function f = mtimes(f, c)
%*	Multiplication of CHEBTECH objects.
%   F*C or C*F multiplies a CHEBTECH F by a scalar or matrix C.
%
%   If F is a verctorised CHEBTECH and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % CHEBTECH * [] = []
    f = []; 
    return
elseif ( ~isa(f, 'chebtech') )       % CHEBTECH is not the first input

    % DOUBLE*CHEBTECH will require that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:CHEBTECH:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a CHEBTECH and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % CHEBTECH * double
    
    % Check dimensions
    if ( size(f.values, 2) ~= size(c, 1) && numel(c) > 1)
        error('CHEBFUN:CHEBTECH:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    if ( numel(c) == 1 )
        f.vscale = f.vscale*abs(c);
        f.epslevel = f.epslevel + eps;
    else
        f.vscale = max(abs(f.values), [], 1);
        % [TODO]: Figure out epslevel.
    end
    
    % If the vertical scale is zero, then set the CHEBTECH to zero:
    if ( all(f.vscale == 0) )
        f.values = zeros(size(f.values, 1), 1);
        f.coeffs = zeros(size(f.values, 1), 1);
    end
    
elseif ( isa(c, 'chebtech') )       % CHEBTECH * CHEBTECH
    
    error('CHEBFUN:CHEBTECH:mtimes:chebtechMtimesFuncheb', ...
        'Use .* to multiply CHEBTECH objects.');
    
else                                % CHEBTECH * ???
    
    error('CHEBFUN:CHEBTECH:mtimes:chebtechMtimesUnknown',...
        'mtimes does not know how to multiply a CHEBTECH and a %s.', class(c));
    
end
