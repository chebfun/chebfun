function f = mtimes(f, c)
%*   Multiplication of CHEBTECH objects.
%   F*C or C*F multiplies a CHEBTECH F by a scalar or matrix C.
%
%   If F is an array-valued CHEBTECH and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % CHEBTECH * [] = []
    f = []; 
    return
elseif ( ~isa(f, 'chebtech') )      % CHEBTECH is not the first input
    % DOUBLE*CHEBTECH requires that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:CHEBTECH:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a CHEBTECH and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % CHEBTECH * double  
    % Check dimensions:
    if ( (size(f.values, 2) ~= size(c, 1)) && (numel(c) > 1) )
        error('CHEBFUN:CHEBTECH:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    if ( numel(c) == 1 )
        % See CHEBTECH CLASSDEF file for documentation on this.
        f.vscale = f.vscale*abs(c);
        f.epslevel = f.epslevel + eps;
    else
        % See CHEBTECH CLASSDEF file for documentation on this.
        f.vscale = max(abs(f.values), [], 1);
        % [TODO]: Figure out vector epslevel.
        f.epslevel = f.epslevel*(1+0*c) + eps;
    end
    
    % If the vertical scale is zero, set the CHEBTECH to zero:
    if ( all(f.vscale == 0) )
        f.values = zeros(1, size(f.values, 2));
        f.coeffs = zeros(1, size(f.values, 2));
    end
    
elseif ( isa(c, 'chebtech') )       % CHEBTECH * CHEBTECH  
    error('CHEBFUN:CHEBTECH:mtimes:chebtechMtimesChebtech', ...
        'Use .* to multiply CHEBTECH objects.');
else                                % CHEBTECH * ???
    error('CHEBFUN:CHEBTECH:mtimes:chebtechMtimesUnknown',...
        'mtimes does not know how to multiply a CHEBTECH and a %s.', class(c));
end

end
