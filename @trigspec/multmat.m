function M = multmat(N, f)
% MULTMAT    Multiplication matrix for the TRIGSPEC class.
%
%  M = MULTMAT(N, F) forms the nxn multiplication matrix
%  representing the multiplication of F in Fourier basis.
% 
%  M = MULTMAT(N, F) also works when F is a vector of Fourier
%  coefficients.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun') ) 
    % Get Fourier coefficients:
    a = get(f, 'coeffs');
elseif ( isa(f, 'bndfun') ) 
    % Get Fourier coefficients:
    a = get(f, 'coeffs');
elseif ( isa( f, 'double') )
    a = f; 
else
    error('TRIGSPEC:ARGIN:TYPE', 'Unrecognised 2nd argument.');
end

% Multiplying by a scalar is easy.
if ( numel(a) == 1 )
    M = a*speye(N);
    return
end

% Deal with even case.
if ( mod(length(a), 2) == 0 )
    a = [a(1)/2; a(2:end); a(1)/2];
end

% Position of constant term.
Na = floor(numel(a)/2) + 1;

% Pad with zeros.
if ( Na < N )
    col = [a(Na:1:end); zeros(N-Na, 1)];
    row = [a(Na:-1:1); zeros(N-Na, 1)];
% Truncate.
else
    col = a(Na:1:end-Na+N);
    row = a(Na:-1:Na-N+1);
end

M = trigspec.sptoeplitz(col, row);

end