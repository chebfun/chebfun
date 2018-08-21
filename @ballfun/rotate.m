function g = rotate(f, lambda, theta)
%ROTATE  Rotate a BALLFUN.
%   ROTATE(f, lambda, theta) is the rotation of f by the azimuthal angle
%   lambda and polar angle theta.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m,n,p] = size(f);

% Evaluation on the Cheb x Fourier x Fourier grid
r = chebpts(m);
lam = pi*trigpts(n);
th = pi*trigpts(p);

G = zeros(m,n,p);
for i = 1:m
    for j = 1:n
        for k = 1:p
            % Transformation to cartesian coordinates
            vec_ijk = sph2cart(lam(j), pi/2-th(k), r(i)).';
            % Rotation of the coordinates
            MrotationZ = rotz(radtodeg(lambda));
            MrotationX = rotx(radtodeg(-theta));
            % Rotation of the normal vector by lambda and theta
            vec_R = MrotationX*MrotationZ*vec_ijk;
            x_0 = vec_R(1); y_0 = vec_R(2); z_0 = vec_R(3); 
            % Transformation to spherical coordinates
            [lam_0, th_0, r_0] = cart2sph(x_0, y_0, z_0);
            th_0 = pi/2-th_0;
            % Evaluation of the function at this point
            G(i,j,k)=feval(f,r_0,lam_0,th_0);
        end
    end
end
% Return the ballfun function
g = ballfun(G, 'vals');
end
