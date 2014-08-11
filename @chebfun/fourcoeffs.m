function varargout = fourcoeffs(f, N)
%FOURCOEFFS   Fourier coefficients of a CHEBFUN.
%   C = FOURCOEFFS(F, N) returns the first N Fourier coefficients of F
%   using complex-exponential form. Specifically: 
%   If N is odd
%       F(x) = C(1)*z^(N-1)/2 + C(2)*z^((N-1)/2-1) + ... + C((N+1)/2) + ... 
%                + C(N)*z^(-(N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(N/2-1) + C(2)*z^(N/2-2) + ... + C(N/2) + ...
%                + C(N-1)*z^(-N/2-1) + 1/2*C(N)*(z^(N/2) + z^(-N/2))
%   where z = exp(1i*omega*x) and omega = 2*pi/L, and L = diff(f.domain). 
%
%   If F is a smooth CHEBFUN (i.e., with no breakpoints), then FOURCOEFFS(F) is
%   equivalent to FOURCOEFFS(F, LENGTH(F)).
%
%   If F is array-valued with M columns, then C is an MxN matrix.
%
%   [A, B] = FOURCOEFFS(F, N) returns the first N Fourier coefficients of F
%   using trignometric form.  Specifically:
%   If N is odd
%      F(x) = A(1)*cos((N-1)/2*omega*x) + B(1)*sin((N-1)/2*omega*x) +  
%             A(2)*cos((N-1)/2-1)*omega*x) + B(2)*sin((N-1)/2-1)*omega*x) + ...
%             + A((N-1)/2)*cos(omega*x) + B((N-1)/2)*sin(omega*x) + A((N+1)/2)
%   If N is even
%      F(x) = A(1)*cos(N/2*omega*x) + B(1)*sin(N/2*omega*x) +  
%             A(2)*cos((N/2-1)*omega*x) + B(2)*sin((N/2-1)*omega*x) + 
%             ... + A(N/2-1)*cos(omega*x) + B(N/2-1)*sin(omega*x) + A(N/2)
%   where omega = 2*pi/L, and L = diff(f.domain). Note that the number of rows
%   in A exceeds the number of rows in B by 1 since A contains the constant
%   term.
%
%   If F is a smooth CHEBFUN (i.e., with no breakpoints), then [A, B] =
%   FOURCOEFFS(F) is equivalent to FOURCOEFFS(F, LENGTH(F)).
%
%   If F is array-valued with M columns, then A and B contain M rows with each
%   row corresponding to the Fourier coefficients for chebfun.
%
% See also CHEBCOEFFS, LEGCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case:
if ( isempty(f) )
    varargout = [];
    return
end

%% Initialize and error checking
numFuns = numel(f.funs);

% If N is not passed in and the numFuns > 1 then throw an error
if ( (nargin == 1) && (numFuns > 1 ) )
    error('CHEBFUN:fourcoeffs:inputN',...
        'Input N is required for piecewise CHEBFUN objects.');
end

if ( nargin == 1 )
    N = length(f);
end

if ( N <= 0 )
    varargout = [];
    return
end
if ( ~isscalar(N) || isnan(N) )
    error('CHEBFUN:fourcoeffs:inputN', 'Input N must be a scalar.');
end

if ( any(isinf(f.domain)) )
    % Fourier coefficients are not allowed for unbounded domains.
    error('CHEBFUN:fourcoeffs:infint', ...
        'Infinite intervals are not supported here.');
end

% Force N to be odd.
N = N + 1 - mod(N,2);

numFuns = numel(f.funs);
if ( numFuns ~= 1 )
    % Attmpet to merge:
    f = merge(f);
    numFuns = numel(f.funs);
end


%% Compute the coefficients.

d = f.domain([1, end]); % Domain of integration
L = diff(d); % Length of the domain

% Modes to compute coefficients.  Need to handle the possible non-symmetry 
% in the modes.
if ( mod(N, 2) == 1 )
    modes = (N-1)/2:-1:-(N-1)/2;
else
    modes = N/2-1:-1:-N/2;
end

if ( numFuns == 1 )
    
    % CHEBCOEFFS() of a smooth piece:
    C = fourcoeffs(f.funs{1}, N).';

else
    % Compute the coefficients via inner products.    
    omega = 2*pi/L;
    x = chebfun('x', d);
    numCols = numColumns(f);
    C = zeros(numCols, N);
    f = mat2cell(f);
    for j = 1:numCols
        coeffsIndex = 1;
        for k = modes
            F = exp(-1i*k*omega*x);
            I = (f{j}.*F);
            C(j, coeffsIndex) = 1/L*sum(I);
            coeffsIndex = coeffsIndex + 1;
        end
    end
end

% The formulas used for computing the coefficients above are explicitly for
% interval [-L/2,L/2] (in the case of one smooth chebfun they are are for
% [-1,1]), not for the domain of f.  We must correct these values to the
% correct interval [a,a+L] by using a change of variables:
%          
% c_n =   \frac{1}{L} \int_{a}^{a+L} f(x) exp(-ikx2\pi/L) dx 
%     =   exp(-ik2\pi(a+L/2)/L) \frac{1}{L} \int_{-L/2}^{L/2} f(t) exp(-ikt2\pi/L) dx 
%
if ( mod(N, 2) == 1 )
    modes = (N-1)/2:-1:-(N-1)/2;
else
    modes = N/2-1:-1:-N/2;
end
change_variables = exp(-1i*modes*2*pi*(d(1) + L/2)/L);
C = bsxfun(@times,C,change_variables);

if ( nargout <= 1 )
    varargout{1} = C;
    
else
    % Return the sign and cosine coefficients for the Fourier series:
    fisOdd = mod(N, 2);
    if ( fisOdd )
        zeroMode = (N+1)/2;
        A = [( C(:,1:zeroMode-1) + C(:,N:-1:zeroMode+1) ), C(:,zeroMode)];
        B = 1i*( C(:,1:zeroMode-1) - C(:,N:-1:zeroMode+1) );
    else  % Non-symmetric case
        zeroMode = N/2;
        A = [C(:,N), ( C(:,1:zeroMode-1) + C(:,N-1:-1:zeroMode+1) ), ...
            C(:,zeroMode)];
        B = 1i*( C(:,1:zeroMode-1) - C(:,N-1:-1:zeroMode+1) );
    end
    if ( isreal(f) )
        A = real(A);
        B = real(B);
    end
    varargout{1} = A;
    varargout{2} = B;
    
end

end

