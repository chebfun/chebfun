function sound(f, varargin)
%SOUND   Play a chebfun as a sound.
%   SOUND(F) overloads the MATLAB sound command for CHEBFUN objects.
%
% See also SING, CHEBTUNE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    error('CHEBFUN:sound:quasi', ...
        'CHEBFUN/SOUND() is not defined for array-valued CHEBFUN objects');
end

if ( isempty(f) )
    f = chebfun(0);
end

[a, b] = domain(f); 

n = 0;
for k = 1:numel(f.funs)
    n = n + length(f.funs{k});
end

Fs = 8192;
n = (b-a)*Fs;
x = linspace(a, b, n);
y = feval(f, x);

if ( numel(varargin) > 0 )
    sound(y, varargin{:});
else
    sound(y, Fs);
end

end
