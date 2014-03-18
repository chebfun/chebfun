function f = diff(f, k, dim)
%DIFF   Derivative of a FOURIERTECH.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the Kth derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this
%   is a finite difference along the columns of an array-valued FOURIERTECH.
%   If F has L columns, an empty FOURIERTECH will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.


%% Check the inputs:

% Trivial case of an empty FOURIERTECH:
if ( isempty(f) )
    return
end

if ( nargin < 2 || isempty(k) )
    % Order of derivative not passed in. Assume 1st derivative by default:
    k = 1; 
elseif ( k == 0 )
    % Nothing to do here!
    return
end    

% Differentiate with respect to the continuous variable by default:
if ( nargin < 3 )
    dim = 1;
end

if ( dim == 1 )
    % Take difference across 1st dimension:
    f = diffContinuousDim(f, k);
else
    % Take difference across 2nd dimension:
    f = diffFiniteDim(f, k);
end

end

function f = diffFiniteDim(f, k)
% Take difference across 2nd dimension.

    % TODO: Tidy and document this.
    if ( k >= size(f, 2) )
        f = f.make();
        return
    else 
        for j = 1:k
            % Differentiate values across dim:
            f.values = diff(f.values, 1, 2);
            % Differentiate coefficients across dim:
            f.coeffs = diff(f.coeffs, 1, 2);
            % Update vscale and epslevel as in PLUS().
            vscale = max(abs(f.values), [], 1);
            ev = f.epslevel.*f.vscale;
            for l = 1:size(f,2)-1
                f.epslevel(l) = ev(l)+ev(l+1);
            end
            f.epslevel(end) = [];
            f.epslevel = f.epslevel./vscale;
            f.vscale = vscale;
        end
    end
end

function f = diffContinuousDim(f, k)
% Differentiate in the first dimension (i.e., df/dx).

    % Get the length:
    N = size(f.coeffs,1);

    % Get the coefficients:
    c = f.coeffs;

    % If n is odd things are easy
    if ( mod(N,2) == 1 )
        % The negative is needed in front of the -(1i)^k term because of
        % the way the coefficients are stored.
        waveNumber = -(-(N-1)/2:(N-1)/2).';
    else
        waveNumber = -([0 (-N/2+1):(N/2-1)]).';
    end
    % Derivative in Fourier space.
    c = bsxfun(@times,c,(1i*waveNumber).^k);
    
    v = f.coeffs2vals(c);

    % [FIXME] Should the epslevel be updated?  This is done in chebfun
    % f.epslevel = n*log(n)*(f.epslevel.*f.vscale);

    % Update the vertical scale
    f.vscale = max(abs(v), [], 1);

    % Store new coefficients and values:
    f.coeffs = c;
    f.values = v;
end