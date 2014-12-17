function f = mtimes(f, c)
%*   Multiplication of TRIGTECH objects.
%   F*C or C*F multiplies a TRIGTECH F by a scalar or matrix C.
%
%   If F is an array-valued TRIGTECH and C is a matrix of appropriate dimension,
%   then the natural matrix multiplication is performed.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(c) )     % TRIGTECH * [] = [].
    
    f = []; 
    return
    
elseif ( ~isa(f, 'trigtech') )      % TRIGTECH is not the first input.
    
    % DOUBLE*TRIGTECH requires that the double is scalar.
    if ( numel(f) > 1 )
        error('CHEBFUN:TRIGTECH:mtimes:size', ...
            'Inner matrix dimensions must agree.');
    end
    
    % C must be a TRIGTECH and F a scalar double. Call MTIMES again.
    f = mtimes(c, f);
    return
    
elseif ( isa(c, 'double') )         % TRIGTECH * double.
    
    % Check dimensions:
    if ( (size(f.values, 2) ~= size(c, 1)) && (numel(c) > 1) )
        error('CHEBFUN:TRIGTECH:mtimes:size2', ...
            'Inner matrix dimensions must agree.');
    end
    
    f.values = f.values*c;
    f.coeffs = f.coeffs*c;
    f.isReal = repmat(all(f.isReal) & isreal(c),1,size(c,2));
    
    if ( numel(c) == 1 )
        % See TRIGTECH CLASSDEF file for documentation on this.
        f.vscale = f.vscale*abs(c);
        f.epslevel = f.epslevel + eps;
    else
        % See TRIGTECH CLASSDEF file for documentation on this.
        vscaleNew = max(abs(f.values), [], 1);
        f.epslevel = ((f.epslevel.*f.vscale)*abs(c))./vscaleNew;
        f.vscale = vscaleNew;

        % Assume condition number 1.
%         glob_acc = max(f.epslevel.*f.vscale);
%         f.vscale = max(abs(f.values), [], 1);
%         f.epslevel = glob_acc./f.vscale;
    end
    
    % If the vertical scale is zero, set the TRIGTECH to zero:
    if ( all(f.vscale == 0) )
        f.values = zeros(1, size(f.values, 2));
        f.coeffs = zeros(1, size(f.values, 2));
    end
    
elseif ( isa(c, 'trigtech') )       % TRIGTECH * TRIGTECH.
    error('CHEBFUN:TRIGTECH:mtimes:trigtechMtimesTrigtech', ...
        'Use .* to multiply TRIGTECH objects.');
    
else                                % TRIGTECH * ???.
    error('CHEBFUN:TRIGTECH:mtimes:trigtechMtimesUnknown',...
        'mtimes does not know how to multiply a TRIGTECH and a %s.', class(c));
end

end