function f = mtimes(f, c)
%*   Multiplication of FOURTECH objects.
%   F*C or C*F multiplies a FOURTECH F by a scalar or matrix C.
%
%   If F is an array-valued FOURTECH and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % FOURTECH * [] = [].
    
    f = []; 
    return
    
elseif ( ~isa(f, 'fourtech') )      % FOURTECH is not the first input.
    
    % DOUBLE*FOURTECH requires that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:FOURTECH:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a FOURTECH and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % FOURTECH * double.
    
    % Check dimensions:
    if ( (size(f.values, 2) ~= size(c, 1)) && (numel(c) > 1) )
        error('CHEBFUN:FOURTECH:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    f.isReal = repmat(all(f.isReal) & isreal(c),1,size(c,2));
    
    if ( numel(c) == 1 )
        % See FOURTECH CLASSDEF file for documentation on this.
        f.vscale = f.vscale*abs(c);
        f.epslevel = f.epslevel + eps;
    else
        % See FOURTECH CLASSDEF file for documentation on this.
        vscaleNew = max(abs(f.values), [], 1);
        f.epslevel = ((f.epslevel.*f.vscale)*abs(c))./vscaleNew;
        f.vscale = vscaleNew;

        % Assume condition number 1.
%         glob_acc = max(f.epslevel.*f.vscale);
%         f.vscale = max(abs(f.values), [], 1);
%         f.epslevel = glob_acc./f.vscale;
    end
    
    % If the vertical scale is zero, set the FOURTECH to zero:
    if ( all(f.vscale == 0) )
        f.values = zeros(1, size(f.values, 2));
        f.coeffs = zeros(1, size(f.values, 2));
    end
    
elseif ( isa(c, 'fourtech') )       % FOURTECH * FOURTECH.
    error('CHEBFUN:FOURTECH:mtimes:fourtechMtimesFourtech', ...
        'Use .* to multiply FOURTECH objects.');
    
else                                % FOURTECH * ???.
    error('CHEBFUN:FOURTECH:mtimes:fourtechMtimesUnknown',...
        'mtimes does not know how to multiply a FOURTECH and a %s.', class(c));
end

end