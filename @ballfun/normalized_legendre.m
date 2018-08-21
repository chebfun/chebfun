function A = normalized_legendre(n)
%NORMALIZED_LEGENDRE   Generate the normalized associated Legendre
% polynomials.
%
% Normalization for the spherical harmonics is \tilde{P}^m_l =
% P^m_l/(2*sqrt(pi)).

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We save the polynomials P^m_l in a matrix.
% In order to create the integral constraints, we have the following
% inequalities
% 0 <= l <= n and 0 <= m <= n
% Then, the discretization for these polynomials is 
% p_tilde = 2*floor(n/2) + 2*p -1 
% If m is the Fourier mode in lambda and l>=m
% then P^m_l will be saved in A(:,l*(l+1)/2+m+1)

% Discretization number of the polynomials
p_tilde = 2*n+1;

% Create the matrix to save the polynomials
A = zeros(p_tilde, n*(n+1)/2+n+1);

% Define useful matrices
trig1 = trigtech( @(t) sin(pi*t) );
Msin = trigspec.multmat(p_tilde, trig1.coeffs );
trig2 = trigtech( @(t) cos(pi*t) );
Mcos = trigspec.multmat(p_tilde, trig2.coeffs );

% Compute the fully normalized associated Legendre polynomials
% P_n^m(cos(theta)), using the Holmes and Featherstone (2002) 
% reference.
for m = 0:n
    if ( m == 0 )
        Pmm = zeros(p_tilde, 1);
        Pmm( floor(p_tilde/2)+1 ) = 1;
    elseif ( m == 1)
        Pmm = sqrt(3)*Msin*Pmm;
    else
        % Compute Pm^m with Pm-1^m-1
        Pmm = sqrt((2*m+1)/(2*m))*Msin*Pmm;
    end
    
    % Save P^m_m
    A(:,m*(m+1)/2+m+1) = (-1)^m*Pmm/sqrt(4*pi);
    
    % Divide by sqrt(2) if m > 0
    if m > 0
        A(:,m*(m+1)/2+m+1) = A(:,m*(m+1)/2+m+1)/sqrt(2);
    end
        
    % Initialize the recurrence (Pm^m-1 does not exist)
    Poldold = zeros(p_tilde, 1);
    Pold = Pmm;

    % Compute P^m_l with the recurrence formula, m+1 <= l <= n
    for l = m+1:n
        anm = sqrt((4*l^2-1)/(l^2-m^2));
        bnm = sqrt((2*l+1)*(l+m-1)*(l-m-1)/((l-m)*(l+m)*(2*l-3)));
        
        % Compute the normalized associated legendre polynomial
        Pl = anm*Mcos*Pold - bnm*Poldold;
        
        % Save it
        A(:, l*(l+1)/2+m+1) = (-1)^m*Pl/sqrt(4*pi);
        
        % Divide by sqrt(2) if m > 0
        if m > 0
            A(:, l*(l+1)/2+m+1) = A(:, l*(l+1)/2+m+1)/sqrt(2);
        end
        
        % Update the polynomials for the recurrence
        Poldold = Pold;
        Pold = Pl;
    end
end
end
