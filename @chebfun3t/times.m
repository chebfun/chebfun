function h = times(f, g, tol)
% .*   CHEBFUN3T multiplication.
%   F.*G multiplies CHEBFUN3 objects F and G. Alternatively F or G could be a
%   double.

if (nargin<3)
    tol = eps;
end

%if ( isa(f, 'chebfun3') )    % CHEBFUN3 .* ???
if ( isa(f, 'chebfun3t') && isa(g, 'chebfun3t'))    % CHEBFUN3 .* ???    
    
%    if ( isa(g, 'double') )  % CHEBFUN2 .* DOUBLE
%        h = mtimes(f, g);
%    elseif ( isa( g, 'chebfun2') )
%        bol = domainCheck(f, g);
%        if ( bol )
            %h = chebfun3t(@(x, y, z) feval(f, x, y, z).*feval(g, x, y, z), f.domain);
            h = chebfun3t(@(x, y, z) feval(f, x, y, z).*feval(g, x, y, z), f.domain, 'eps', tol);
%        else
%            error('CHEBFUN:CHEBFUN2:times:domain', 'Inconsistent domains');
%        end
    else
        error('CHEBFUN:CHEBFUN3T:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
    
%else
%    h = times(g, f);
%end

end

function hCoeffs = coeff_times(fCoeffs, gCoeffs)
%COEFF_TIMES(FC, GC)   Multiplication in coefficient space
%   HC = COEFF_TIMES(FC, GC) returns the vector of Chebyshev coefficients, HC,
%   resulting from the multiplication of two functions with FC and GC
%   coefficients. The vectors have already been prolonged.

%   Multiplication in coefficient space is a Toeplitz-plus-Hankel-plus-rank-one
%   operator (see Olver & Townsend, A fast and well-conditioned spectral method,
%   SIAM Review, 2013). This can be embedded into a Circular matrix and applied
%   using the FFT:

mn = length(fCoeffs);
t = [2*fCoeffs(1,:) ; fCoeffs(2:end,:)];                    % Toeplitz vector.
x = [2*gCoeffs(1,:) ; gCoeffs(2:end,:)];                    % Embed in Circulant.
xprime = fft([x ; x(end:-1:2,:)]);                % FFT for Circulant mult.
aprime = fft([t ; t(end:-1:2,:)]);
Tfg = ifft(aprime.*xprime);                   % Diag in function space.
hCoeffs = .25*[Tfg(1,:); Tfg(2:end,:) + Tfg(end:-1:2,:)];% Extract out result.
hCoeffs = hCoeffs(1:mn,:);                                % Take the first half.
end

function hCoeffs = coeff_times3D(fCoeffs, gCoeffs)
[m_f, n_f, p_f] = size(fCoeffs);
[m_g, n_g, p_g] = size(gCoeffs);
hCoeffs = zeros(m_f+m_g-1, n_f+n_g-1, p_f+p_g-1);

hVals = chebfun3t.coeffs2vals(fCoeffs, gCoeffs);
hCoeffs = chebfun3t.vlas2coeffs(hVals);

mn = length(fCoeffs);
t = [2*fCoeffs(1, :); fCoeffs(2:end, :)];                    % Toeplitz vector.
x = [2*gCoeffs(1, :); gCoeffs(2:end, :)];                    % Embed in Circulant.
xprime = fft([x; x(end:-1:2, :)]);                % FFT for Circulant mult.
aprime = fft([t; t(end:-1:2, :)]);
Tfg = ifft(aprime .* xprime);                   % Diag in function space.
hCoeffs = .25*[Tfg(1, :); Tfg(2:end, :) + Tfg(end:-1:2, :)];% Extract out result.
hCoeffs = hCoeffs(1:mn, :);                                % Take the first half.
end