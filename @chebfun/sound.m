function sound(f, varargin)
%SOUND   Play a CHEBFUN as a sound.
%   SOUND(F) overloads the MATLAB SOUND command for CHEBFUN objects.
%
% See also SING, CHEBTUNE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document this?

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:sound:quasi', ...
        'CHEBFUN/SOUND() is not defined for array-valued CHEBFUN objects');
end

if ( isempty(f) )
    f = chebfun(0);
end

[a, b] = domain(f); 

Fs = 8192;
n = (b - a)*Fs;
x = linspace(a, b, n);
y = feval(f, x);

if ( numel(varargin) > 0 )
    sound(y, varargin{:});
else
    sound(y, Fs);
end

end
