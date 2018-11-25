% function f = randnfunball(lambda)
% %RANDNFUNBALL   Smooth random function on the unit sphere
% %   F = RANDNFUNBALL(LAMBDA) returns a SPHEREFUN of maximum
% %   wavelength of about 2pi/LAMBDA and standard normal distribution
% %   N(0,1) at each point.  F is obtained from a combination of
% %   spherical harmonics with random coefficients.
% %
% %
% % See also RANDNFUN3, RANDNFUNSPHERE, RANDNFUNDISK.
% 
% % Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% % See http://www.chebfun.org/ for Chebfun information.
% 
% % Parse the inputs
% if nargin == 0
%     lambda = 1;
% end
% 
% % Since the unit sphere has circumference L = 2*pi, the following
% % choice matches the choice deg = floor(L/lambda) in randnfun.
% deg = floor(2*pi/lambda);
% 
% % The coefficient for the SH of degree l and order m is stored in 
% % c(l^2+l+m+1)
% c = randn((deg+1)^2,1);
% c = sqrt(4*pi/nnz(c))*c;         % normalize so variance is 1
% 
% % TODO : remove this
% %deg = 5;
% %c = zeros(6^2);
% %c(5^2+5-3+1) = 1;
% 
% f = ballfun(sumSolHarm(deg, c), 'coeffs');
% end
% 
% function F = sumSolHarm(deg, c)
% %NORMALIZED_LEGENDRE   Generate the fully normalized associated Legendre
% % polynomials P^m_max_l_max
% 
% % Copyright 2018 by The University of Oxford and The Chebfun Developers.
% % See http://www.chebfun.org/ for Chebfun information.
% 
% % Initialize F in order r, th, lambda
% S = [deg+1,2*deg+1,2*deg+1];
% F = zeros(S);
% 
% % Discretization number of the polynomial
% p = 2*deg+1;
% 
% % Interpolation points
% th = pi*trigpts(p);
% 
% % Precompute vector cos(th), sin(th)^m and r^l
% r = chebtech2.vals2coeffs(chebpts(deg+1, [-1 1]).^(0:deg));
% CosTh = cos(th);
% SinTh = sin(th).^(0:deg);
% 
% % Initialize P^0_0 / u^0
% Pmm = ones(p, 1);
% Pinitm = ones(p, 1);
% 
% %% Compute P^m_m/u^m
% for m = 0:deg
%     % Compute Pm^m with Pm-1^m-1
%     Pmm = sqrt((2*m+1)/(2*m-(m==1)+(m==0)))*Pmm;
%     
%     %% Build the SH
%     Plm = (-1)^m*SinTh(:,m+1).*Pmm/sqrt(4*pi*(1+(m>0)));
%     % FFT to recover the coefficients
%     Plm = trigtech.vals2coeffs(Plm);
%     F(:,:,deg+m+1) = F(:,:,deg+m+1) + c((m+1)^2)*build_ballfun_Legendre(r(:,m+1), Plm, m, m);
%     if m ~= 0
%         F(:,:,deg-m+1) = F(:,:,deg-m+1) + c(m^2+1)*build_ballfun_Legendre(r(:,m+1), Plm, m, -m);
%     end
%     
%     Poldold = Pinitm;
%     Pold = Pmm;
%     
%     %% Compute P^m_l / u^m with the recurrence formula, m_max+1 <= l <= l_max
%     for l = m+1:deg
%         anm = sqrt((4*l^2-1)/((l-m)*(l+m)));
%         bnm = sqrt((2*l+1)*(l+m-1)*(l-m-1)/((l-m)*(l+m)*(2*l-3)));
%         % Compute the normalized associated legendre polynomial P^m_l/u^m
%         Pl = anm*CosTh.*Pold - bnm*Poldold;
%         
%         %% Build the SH
%         Plm = (-1)^m*SinTh(:,m+1).*Pl/sqrt(4*pi*(1+(m>0)));
%         % FFT to recover the coefficients
%         Plm = trigtech.vals2coeffs(Plm);
%         F(:,:,deg+m+1) = F(:,:,deg+m+1) + c(l^2+l+m+1)*build_ballfun_Legendre(r(:,l+1), Plm, l, m);
%         if m ~= 0
%             F(:,:,deg-m+1) = F(:,:,deg-m+1) + c(l^2+l-m+1)*build_ballfun_Legendre(r(:,l+1), Plm, l, -m);     
%         end
%         
%         %% Update the polynomials for the recurrence
%         Poldold = Pold;
%         Pold = Pl;
%     end
%     
% end
% 
% % Permute F to get r, lambda, theta
% F = permute(F, [1,3,2]);
% 
% end
% 
% function F = build_ballfun_Legendre(monomial, Plm, l, m)
% % Take the coeffs of P^|m|_l and the coeffs of r^l and return the tensor of 
% % coeffs of r^lP^m_l 
% 
% if m < 0
%    Plm = (-1)^m*Plm; 
% end
% 
% % Normalize the solid harmonic so that its two-norm over the ball is 1
% % This is the normalization constant so that the integral of r^l is 1
% Plm = Plm*sqrt(2*l+3);
% 
% % Return the Spherical Harmonic Y^m_l
% F = monomial*Plm.';
% end

%% Old version
function f = randnfunball(lambda)
%RANDNFUNBALL   Smooth random function on the unit ball
%   F = RANDNFUNBALL(LAMBDA) returns a BALLFUN of maximum
%   frequency about 2pi/LAMBDA and standard normal distribution N(0,1)
%   at each point.  F is obtained by restricting output of RANDNFUN3
%   to the unit ball.
%
%   RANDNFUNBALL() uses the default value LAMBDA = 1.
%
% See also RANDNFUN3, RANDNFUNSPHERE, RANDNFUNDISK.

% Copyright 2018 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 0
    lambda = 1;
end

% Create a randnfun3 function
fcart = randnfun3( lambda );

% Create a tensor product grid in (r,lam,th)
sz = min(length(fcart),50);

% Sample on a [sz,sz,sz] grid:
r = chebpts( sz );
lam  = pi*trigpts( 2*sz+1 );
th  = pi*trigpts( 2*sz+1 );

[rr, ll, tt] = ndgrid(r, lam, th);

% Sample the random function over the square and create the ballfun:
xx = rr.*sin(tt).*cos(ll);
yy = rr.*sin(tt).*sin(ll);
zz = rr.*cos(tt);
F = feval(fcart, xx, yy, zz);

% Simplify and return the output
f = simplify(ballfun(F));

end