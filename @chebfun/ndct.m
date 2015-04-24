function v = ndct( c )
%NDCT  Evaluate a Chebyshev series at Gauss-Legendre nodes.
%
%   V = NDCT(C) returns a column vector V such that
%       V(k) = C(1)*T_0(x(k)) + C(2)*T_1(x(k)) + ... C(N)*T_{N-1}(x(k)), 
%   where T_j(x) is the degree j Chebyshev polynomial and x is the vector of
%   Gauss-Legendre nodes (as returned by x = legpts(size(C,1))).
%
% See also CHEBFUN.DLT, CHEBFUN.IDLT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: Uses the 'cheb_1' version of the algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = eps;                                    % Tolerance
L = 18;                                       % Max number of terms in series

[n, m] = size(c); nn = (0:n-1).';             % Transform size
[~, ~, ~, t_leg] = legpts(n);                 % Legendre grid (in theta)
t_leg = t_leg(end:-1:1);                      % Get right-left ordering
t_cheb = ((0:n-1).'+.5)*pi./n;                % Chebyshev nodes (in theta space)
dt = (t_leg - t_cheb);                        % (t_leg - t_cheb)
dt(n:-1:ceil((n+1)/2)) = -dt(1:ceil(n/2));    % More accurate near t = 0
dt = repmat(dt, 1, m); nn = repmat(nn, 1, m); % Support matrix input.
Dt = 1; scl = 1;                              % Initialise
v = dct3_scaled(c);                           % ell = 0 term
for l = 1:(L-1)                               % Remaining terms in sum
    c = nn.*c; Dt = Dt.*dt; scl = scl*l;      % Update terms
    sgn = (-1)^floor((l+1)/2);                % +--++--
    if ( mod(l, 2) )
        dv = (sgn/scl)*Dt.*dst3_shifted(c);   % Sine terms
    else       
        dv = (sgn/scl)*Dt.*dct3_scaled(c);    % Cosine terms
    end
    v = v + dv;
    if ( norm(dv, inf) < tol ), break, end
end
v = v(end:-1:1,:);                            % Return to left-right ordering

% Ensure real/imag output for real/imag input:
if ( isreal(c) ), v = real(v); elseif ( isreal(1i*c) ), v = imag(v); end

end

function v = dct3_scaled(c)
% Scaled DCT-III.
c(1,:) = 2*c(1,:);                            % Scale first coefficient
v = chebfun.dct(c, 3);                        % Compute DCT-III
end

function v = dst3_shifted(c)
% Shifted DST-III.
c = [c(2:end,:) ; zeros(1, size(c, 2))];      % Shift first coefficient
v = chebfun.dst(c, 3);                        % Compute DST-III.
end
