function f = times(f, g, varargin)
%.*   CHEBTECH multiplication.
%   F.*G multiplies CHEBTECH objects F and G or a CHEBTECH by a scalar if either
%   F or G is a scalar.
%
%   If F is an array-valued CHEBTECH, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MTIMES, RDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% CHEBTECH * [] = []:
if ( isempty(f) || isempty(g) )
    f = [];
    return
end

if ( ~isa(f, 'chebtech') )      % Ensure F is a CHEBTECH
    f = times(g, f, varargin{:});
    return
elseif ( isa(g, 'double') )     % CHEBTECH .* double
    
    % Do the multiplication:
    if ( size(g, 2) > 1 )
        f.coeffs = bsxfun(@times, f.coeffs, g);
    else
        f.coeffs = f.coeffs*g;
    end
    f.epslevel = f.epslevel + eps(g);
    return
    
elseif ( size(f.coeffs, 1) == 1 )
    % If we have (constant CHEBTECH).*CHEBTECH, reverse the order and call TIMES
    % again:
    f = times(g, f.coeffs);
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
elseif ( size(g.coeffs, 1) == 1)
    % If we have CHEBTECH.*(constant CHEBTECH), convert the (constant CHEBTECH)
    % to a scalar and call TIMES again:
    f = times(f, g.coeffs);
    f.epslevel = max(f.epslevel, g.epslevel);
    return
end

% Get the size of each CHEBTECH:
[fn, fm] = size(f.coeffs);
[gn, gm] = size(g.coeffs);

fNew = flipud(f.coeffs);
gNew = flipud(g.coeffs);

% prolong:
fNew((fn+1):(fn+gn+1),:) = 0;
gNew((gn+1):(fn+gn+1),:) = 0;

% Check dimensions:
if ( fm ~= gm )
    if ( fm == 1 )
        % Allow [Inf x 1] .* [Inf x m].
        fNew = repmat(fNew, 1, gm);
    elseif ( gm == 1 )
        % Allow [Inf x m] .* [Inf x 1].
        gNew = repmat(gNew, 1, fm);
    else
        error('CHEBFUN:CHEBTECH:times:dim2', ...
            'Inner matrix dimensions must agree.');
    end
end

% Check for two cases where the output is known in advance to be positive,
% namely F == conj(G) or F == G and isreal(F).
pos = false;

% Multiply values:
if ( isequal(f, g) )
    coeffs = coeff_times( fNew, gNew );
    if ( isreal(f) )
        pos = true;
    end
elseif ( isequal(conj(f), g) )
    coeffs = coeff_times( conj(fNew), gNew );
    pos = true;
else
    coeffs = coeff_times( fNew, gNew );
end

% Assign values and coefficients back to f:
coeffs = flipud(coeffs);

f.coeffs = coeffs; 

% Update vscale, epslevel, and ishappy:
vscale = getvscl(f);
% See CHEBTECH CLASSDEF file for documentation on this:
f.epslevel = (f.epslevel + g.epslevel) .* (vscale.*g.vscale./vscale);
% f.vscale  = vscale;
f.ishappy = f.ishappy && g.ishappy;

% Simplify!
f = simplify(f);

if ( pos )
    % Here we know that the product of F and G should be positive. However,
    % SIMPLIFY may have destroyed this property, so we enforce it.
    values = f.coeffs2vals(f.coeffs); 
    values = abs(values);
    f.coeffs = f.vals2coeffs(values);
end

end


function hc = coeff_times(fc, gc)
% COEFF_TIMES(FC, GC) multiplication in coefficient space
% 
% HC = COEFF_TIMES(FC, GC) returns the vector of Chebyshev coefficients, HC,
% resulting from the multiplication of two functions with FC and GC
% coefficients. The vectors have already been prolonged. 

% Multiplication in coefficient space is a Toeplitz-plus-Hankel-plus-rank-one
% operator (see Olver & Townsend, A fast and well-conditioned spectral
% method, SIAM Review, 2013). This can be embedded into a Circular matrix and
% applied using the FFT: 

mn = length(fc);
t = [2*fc(1,:) ; fc(2:end,:)];                    % Toeplitz vector.
x = [2*gc(1,:) ; gc(2:end,:)];                    % Embed in Circulant.
xprime = fft([x ; x(end:-1:2,:)]);              % FFT for Circulant mult.
aprime = fft([t ; t(end:-1:2,:)]);
Tfg = ifft(aprime.*xprime);                   % Diag in function space.
hc = .25*[Tfg(1,:); Tfg(2:end,:) + Tfg(end:-1:2,:)];% Extract out result.
hc = hc(1:mn,:);                                % Take the first half.

end
