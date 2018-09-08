function f = solharm(l,m)
% SOLHARM Complex-valued, solid harmonic of degree L, order M.
%
%   Y = SPHHARM(L, M) returns the degree L, order M real-valued 
%   solid harmonic on the ball.  Y is normalized so that its two-norm
%   over the sphere is 1.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( l < abs(m) )
    error('CHEBFUN:BALLFUN:solHarm', 'Degree must be >= order for solid harmonic ');
else
    % abs(m) <= l
    % Compute P^abs(m)_l 
    S = [l+1,2*abs(m)+1,2*l+1];
    Plm = normalized_legendre(l,abs(m));
    if m < 0
       Plm = (-1)^m*Plm; 
    end
    % Normalize the solid harmonic so that its two-norm over the ball is 1
    % This is the normalization constant so that the integral of r^l is 1
    Plm = Plm*sqrt(2*l+3);
    
    % Compute the chebyshev coefficients of r^l
    r = chebpts(l+1);
    monomial = r.^l;
    monomial = chebtech2.vals2coeffs(monomial);
    
    % Return the Spherical Harmonic Y^m_l
    F = zeros(S);
    F(:,abs(m)+m+1,:) = reshape(monomial*Plm.',l+1,1,2*l+1);
    f = ballfun(F,'coeffs');
end

end


function Pold = normalized_legendre(l_max,m_max)
%NORMALIZED_LEGENDRE   Generate the fully normalized associated Legendre
% polynomials P^m_max_l_max
%
% Normalization for the spherical harmonics is \tilde{P}^m_l =
% P^m_l/(2*sqrt(pi)).

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Below is the implementation the Modified Forward Column (MFC) method 
% described in the Holmes and Featherstones paper (2002)
% It computes P^m_n/u^m by a stable recurrence to avoid numerical errors
% near the poles

% Discretization number of the polynomial
p = 2*l_max+1;

% Interpolation points
th = pi*trigpts(p);
% Precompute the cos(th)
CosTh = cos(th);

%% Compute P_m_max^m_max / u^m_max
for m = 0:m_max
    if ( m == 0 )
        Pmm = ones(p, 1);
    elseif ( m == 1)
        Pmm = sqrt(3)*Pmm;
    else
        % Compute Pm^m with Pm-1^m-1
        Pmm = sqrt((2*m+1)/(2*m))*Pmm;
    end
end

% Initialize the recurrence (Pm^m-1 does not exist)
Poldold = zeros(p, 1);
Pold = Pmm;

%% Compute P^m_l with the recurrence formula, m_max+1 <= l <= l_max
for l = m_max+1:l_max
    anm = sqrt((4*l^2-1)/((l-m_max)*(l+m_max)));
    bnm = sqrt((2*l+1)*(l+m_max-1)*(l-m_max-1)/((l-m_max)*(l+m_max)*(2*l-3)));
    % Compute the normalized associated legendre polynomial
    Pl = anm*CosTh.*Pold - bnm*Poldold;

    % Update the polynomials for the recurrence
    Poldold = Pold;
    Pold = Pl;
end

% Normalize the polynomial and recover associated Legendre polynomials
% Divide by sqrt(2) if m_max > 0
Pold = (-1)^m*sin(th).^m.*Pold/sqrt(4*pi*(1+(m>0)));

% FFT to recover the coefficients
Pold = trigtech.vals2coeffs(Pold);
end