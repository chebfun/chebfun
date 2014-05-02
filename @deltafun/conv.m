function h = conv(f, g)
%CONV   Convolution of DELTAFUN objects.
%   H = CONV(F, G) produces the convolution of DELTAFUN objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                    /
%                   -
%   Example:
%     f = chebfun(1/2); g = f;
%     subplot(2, 2, 1), plot(f)
%     for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%     figure, for j = 1:4, subplot(2,2,j), plot(g), g = diff(g); end

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = deltafun();
    return
end

if ( ~isa(f, 'deltafun') )
    % Then g must be a deltafun
    h = conv(g, f);
    return
end

h = conv(f.funPart, g.funPart);

warning('CHEBFUN:deltafun:noSupport', ...
    'DELTAFUN doesn''t fully support CONV() yet.');
  
end

