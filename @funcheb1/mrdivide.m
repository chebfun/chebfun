function f = mrdivide(f, c)
%/	Right scalar divide for a FUNCHEB1.
%   F/C divides the FUNCHEB1 F by a scalar C.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO] Add support for vectorised FUNCHEB1 objects.
% [TODO] What if c == 0 ?

if ( isa(c,'double') )
    
    if ( numel(c) > 1 )
        error('CHEBFUN:FUNCHEB1:mrdivide:size', ...
            'No support for FUNCHEB1 mrdivide with vectors or matrices.');
    end
    
    f.values = f.values/c;      % Divide values
    f.coeffs = f.coeffs/c;      % Divide coeffs
    f.vscale = f.vscale/norm(c,inf); % Divide vscale
    
else
    
    error('FUNCHEB1:mrdivide:funfun', 'Use ./ to divide by a FUNCHEB1.');
    
end

end