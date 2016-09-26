function I = sum3(f)
%SUM3   Triple definite integral of a CHEBFUN3T over [a, b] x [c, d] x [e g].
%   This is a direct generalization of Clenshaw-Curtis quadrature as 
%   implemented in Page 77 of Battles' PhD thesis.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) ) 
    I = []; 
    return; 
end

coeffs = f.coeffs;
% Rescaling factors:
rescaleFactor = 0.5*diff(f.domain(1:2));                 % (b - a)/2
rescaleFactor = rescaleFactor * 0.5*diff(f.domain(3:4)); % (d - c)/2
rescaleFactor = rescaleFactor * 0.5*diff(f.domain(5:6)); % (g - e)/2

if ( nargin == 1 )
    [nx, ny, nz] = size(coeffs);
    I1 = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            I1(i,j) = squeeze(coeffs(i,j,1:2:nz)).'*(2./(1-(0:2:nz-1)'.^2));
        end
    end
    I2 = zeros(nx,1);
    for i = 1:nx
        I2(i) = I1(i,1:2:ny)*(2./(1-(0:2:ny-1)'.^2));
    end
    I = I2(1:2:nx).'*(2./(1-(0:2:nx-1)'.^2));
end

% Rescale the output:
I = I*rescaleFactor;

end