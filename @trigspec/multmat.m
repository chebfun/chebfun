function M = multmat(N, f)
% MULTMAT    Multiplication matrix for the TRIGSPEC class.
%
%  M = MULTMAT(N, F) forms the nxn multiplication matrix
%  representing the multiplication of F in Fourier basis.
% 
%  M = MULTMAT(N, F) also works when F is a vector of Fourier
%  coefficients.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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

% Coefficient of index k=0.
Na = floor(numel(a)/2);

% Pad with zeros.
if ( Na < N )
    row = [a(Na+1:-1:1).', zeros(1, N-Na-1)];
    column = [a(Na+1:1:end); zeros(N-Na-1, 1)];
% Truncate.
else
    row = a(1+Na:-1:Na-N+2).';
    column = a(1+Na:1:end-Na+N-1);
end
M = trigspec.sptoeplitz(row, column);

end