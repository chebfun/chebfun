function f = mtimes(f, c)
%*   Multiplication of FOURIERTECH objects.
%   F*C or C*F multiplies a FOURIERTECH F by a scalar or matrix C.
%
%   If F is an array-valued FOURIERTECH and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % FOURIERTECH * [] = []
    f = []; 
    return
elseif ( ~isa(f, 'fouriertech') )      % FOURIERTECH is not the first input
    % DOUBLE*FOURIERTECH requires that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:FOURIERTECH:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a FOURIERTECH and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % FOURIERTECH * double  
    % Check dimensions:
    if ( (size(f.values, 2) ~= size(c, 1)) && (numel(c) > 1) )
        error('CHEBFUN:FOURIERTECH:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    if ( numel(c) == 1 )
        % See FOURIERTECH CLASSDEF file for documentation on this.
        f.vscale = f.vscale*abs(c);
        f.epslevel = f.epslevel + eps;
    else
        % See FOURIERTECH CLASSDEF file for documentation on this.
        vscaleNew = max(abs(f.values), [], 1);
        f.epslevel = ((f.epslevel.*f.vscale)*abs(c))./vscaleNew;
        f.vscale = vscaleNew;

        % Assume condition number 1.
%         glob_acc = max(f.epslevel.*f.vscale);
%         f.vscale = max(abs(f.values), [], 1);
%         f.epslevel = glob_acc./f.vscale;
    end
    
    % If the vertical scale is zero, set the FOURIERTECH to zero:
    if ( all(f.vscale == 0) )
        f.values = zeros(1, size(f.values, 2));
        f.coeffs = zeros(1, size(f.values, 2));
    end
    
elseif ( isa(c, 'chebtech') )       % FOURIERTECH * FOURIERTECH  
    error('CHEBFUN:FOURIERTECH:mtimes:chebtechMtimesChebtech', ...
        'Use .* to multiply FOURIERTECH objects.');
else                                % FOURIERTECH * ???
    error('CHEBFUN:FOURIERTECH:mtimes:chebtechMtimesUnknown',...
        'mtimes does not know how to multiply a FOURIERTECH and a %s.', class(c));
end

end
