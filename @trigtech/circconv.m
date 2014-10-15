function f = circconv(f, g)
%CONV   Circular convolution of TRIGTECH objects.
%   H = CIRCCONV(F, G) produces the convolution of TRIGTECH objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [-pi, pi]
%                    /
%                   -
%   Note that CIRCCONV only supports smooth periodic functions on [-pi,pi].
%
%   Example:
%     f = trigtech(@(x) exp(cos(40*pi*x))); 
%     g = trigtech(@(x) exp(-(20*x).^2);
%     h = circconv(f,g);
%     plot(h);

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    f = trigtech();
    return
end

% No support for array-valued trigtech objects:
if ( (size(f, 2) > 1) || (size(g, 2) > 1) )
    error('CHEBFUN:TRIGTECH:conv:array', ...
        'No support for array-valued TRIGTECH objects.');
end

% Get the sizes of the TRIGTECH objects
nf = size(f.coeffs, 1);
ng = size(g.coeffs, 1);

% Make the TRIGTECH objects the same length.
if ( nf > ng )
    % Increase the length of g (via PROLONG):
    g = prolong(g, nf);
elseif ( nf < ng )
    % Increase the length of f (via PROLONG):
    f = prolong(f, ng);
end
n = size(f.coeffs,1);

% Convolution is just multiplication of the Fourier coefficients.
% Shift g horizontally to -1.
g = circshift(g,-1);
f.values = 2/n*ifft(fft(f.values).*fft(g.values));
f.coeffs = f.vals2coeffs(f.values);

% Update the vscale.
f.vscale = max(abs(f.values), [], 1);

% Scale the epslevel relative to the largest column:
vscale = f.vscale;
f.epslevel = 10*eps(max(f.vscale));
vscale(vscale <= f.epslevel) = 1;
f.epslevel = f.epslevel./vscale;

f = simplify(f);

f.ishappy = f.ishappy && g.ishappy;
f.isReal = f.isReal && g.isReal;  % Are you real happy though?

f.values(:,f.isReal) = real(f.values(:,f.isReal));

if ( f.ishappy )
    f = simplify(f, f.epslevel);
end

end