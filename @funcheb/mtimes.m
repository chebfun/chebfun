function f = mtimes(f, c)
%*	Multiplication of FUNCHEB objects.
%   F*C or C*F multiplies a FUNCHEB F by a scalar or matrix C.
%
%   If F is a verctorised FUNCHEB and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % FUNCHEB * [] = []
    f = []; 
    return
elseif ( ~isa(f, 'funcheb') )       % FUNCHEB is not the first input

    % DOUBLE*FUNCHEB will require that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:FUNCHEB:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a FUNCHEB and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % FUNCHEB * double
    
    % Check dimensions
    if ( size(f.values, 2) ~= size(c, 1) && numel(c) > 1)
        error('CHEBFUN:FUNCHEB:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    f.vscale = f.vscale*abs(c);
    
    % If the vertical scale is zero, then set the FUNCHEB to zero:
    if ( all(f.vscale == 0) )
        f.values = zeros(size(f.values, 1), 1);
        f.coeffs = zeros(size(f.values, 1), 1);
    end
    
elseif ( isa(c, 'funcheb') )       % FUNCHEB * FUNCHEB
    
    error('CHEBFUN:FUNCHEB:mtimes:funchebMtimesFuncheb', ...
        'Use .* to multiply FUNCHEB objects.');
    
else                                % FUNCHEB * ???
    
    error('CHEBFUN:FUNCHEB:mtimes:funchebMtimesUnknown',...
        'mtimes does not know how to multiply a FUNCHEB and a %s.', class(c));
    
end
