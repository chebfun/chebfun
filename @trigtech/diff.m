function f = diff(f, k, dim)
%DIFF   Derivative of a TRIGTECH.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the Kth derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this
%   is a finite difference along the columns of an array-valued TRIGTECH.
%   If F has L columns, an empty TRIGTECH will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% Trivial case of an empty TRIGTECH:
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
% Take kth difference across 2nd dimension (i.e., across columns).

    if ( k >= size(f, 2) )
        % The output will be an empty TRIGTECH:
        f = f.make();
    else 
        for j = 1:k
            % Differentiate values across dim:
            f.values = diff(f.values, 1, 2);
            % Differentiate coefficients across dim:
            f.coeffs = diff(f.coeffs, 1, 2);
            f.isReal = ~logical(diff(~f.isReal));
        end
    end
end

function f = diffContinuousDim(f, k)
% Differentiate in the first dimension (i.e., df/dx).

    % Get the length:
    N = size(f.coeffs,1);

    % Get the coefficients:
    c = f.coeffs;
    
    % Code depends on parity of N
    if ( mod(N,2) == 1 ) % Odd
        % Since N is odd things are easy.
        waveNumber = (-(N-1)/2:(N-1)/2).';
    else
        % Wave numbers are non-symmetric in the even case
        waveNumber = (-N/2:N/2-1).';
    end 

    % Derivative in Fourier space.
    c = bsxfun(@times,c,(1i*pi*waveNumber).^k);
    
    v = f.coeffs2vals(c);
    v(:,f.isReal) = real(v(:,f.isReal));

    % Store new coefficients and values:
    f.coeffs = c;
    f.values = v;
end
