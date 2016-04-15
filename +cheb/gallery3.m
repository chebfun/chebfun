function varargout = gallery3(name)
%CHEB.GALLERY3   CHEBFUN3 example functions.
%
%   [F, FA] = CHEB.GALLERY3(NAME) also returns the anonymous function FA 
%   used to define the function. 
%   CHEB.GALLERY3 with no output argument creates a plot of the selected
%   function:
%
%   genz1        An oscillatory integration test function
%   genz2        An integration test function with a corner peak
%   genz3        Another integration test function with a Gaussian peak
%   genz4        An integration test function with a product peak
%   griewank     A global minimization problem with many widespread local 
%                minima. It's global minimum is at origin with value zero.
%  rastrigin     Another global minimization test problem with lots of
%                local minima. It's global minimum is at origin with value zero.
%   runge        A 3D Runge function
%   wagon        A challenging optimization problem suggested by Wagon
%   shubert      Shubert function from global optimization
%   rosenbrock   Rosenbrock function from global optimization. It's global
%                minimum is zero at (1,1,1).
%   levy         Levy function. Another global minimization test problem
%                with minimum value zero at (1,1,1).
%   aurentz      A 3D radial bump function suggested by Aurentz
%   clebsch      Clebsch Diagonal Cubic from algebraic geometry
%   cassini      Cassini Surface from algebraic geometry
%   barth10      Barth Decic from algebraic geometry
%   wave1        A standing wave function
%   wave2        A travelling wave function
%   slevinsky    A function related to the heat and wave equations 
%                suggested by Slevinsky
%   Also, the following three functions 
%   3dxy, 
%   4px,
%   32m1,
%   5gz3y, 
%   are a few of the complete wave functions of the Hydrogen atom which 
%   represent its orbitals.
%
%   Gallery functions are subject to change in future releases of Chebfun.
%
% See also CHEB.GALLERY, CHEB.GALLERYTRIG and CHEB.GALLERY2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If the user did not supply an input, return a random function from the
% gallery.
if ( nargin == 0 )
    names = {'wagon', 'runge', 'clebsch', 'cassini', 'barth10', ...
        'shubert', 'rosenbrock', 'aurentz', 'genz1', 'genz2', 'genz3', ...
        'genz4', 'griewank', 'rastrigin', 'wave1', 'wave2', 'slevinsky', ...
        '3dxy', '4px', '32m1', '5gz3y'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)

   case 'genz1'
       % http://www.sfu.ca/~ssurjano/oscil.html
       fa = @(x,y,z) cos(pi + 5 * (x+y+z));
       f = chebfun3(fa);
       % Try e.g., plot(sum(f))
       
   case 'genz2'
       %http://www.sfu.ca/~ssurjano/copeak.html
       fa = @(x,y,z) 1./(16 +5 * (x+y+z));
       f = chebfun3(fa);
       % Try e.g., surf(f)
       
   case 'genz3'
       % http://www.sfu.ca/~ssurjano/gaussian.html
       fa = @(x,y,z) exp(-25*(x.^2+y.^2+z.^2));
       f = chebfun3(fa);
       % Try e.g., surf(f)

    case 'genz4'
       % http://www.sfu.ca/~ssurjano/prpeak.html
       fa = @(x,y,z) 1./((0.04+x.^2) .* (0.04+y.^2) .* (0.04+z.^2));
       f = chebfun3(fa);
       % Try e.g., surf(f)

       case 'griewank'
       %http://www.sfu.ca/~ssurjano/griewank.html
       fa = @(x,y,z) x.^2/4000 + y.^2/4000 + z.^2/4000 - ...
           cos(x).*cos(y/sqrt(2)).*cos(z/sqrt(3)) + 1;
       f = chebfun3(fa, [-50 50 -50 50 -50 50]);
       % Try e.g., plot(f) or surf(f)    
    
    case 'wagon'
        % Function from Page 99 of F. Bornemann, D. Laurie, S. Wagon and J. 
        % Waldvogel, The SIAM 100-Digit Challenge, SIAM, 2004.
        fa = @(x,y,z) exp(sin(50*x)) + sin(60*exp(y)).*sin(60*z) + ...
            sin(70*sin(x)).*cos(10*z) + sin(sin(80*y)) - sin(10*(x+z)) +...
            (x.^2 + y.^2 + z.^2)/4;
        f = chebfun3(fa);
        
    case 'runge'
        % A 3D Runge function:        
        fa = @(x,y,z) 1./(1+x.^2+y.^2+z.^2);
        f = chebfun3(fa);
        
    case 'shubert'
        %Shubert Function 
        % http://www.sfu.ca/~ssurjano/shubert.html
        fa = @(x,y,z) (cos(2*x+1) + 2*cos(3*x+2) + 3*cos(4*x+3) + ...
            4*cos(5*x+4) + 5*cos(4*x+5)) .* (cos(2*y+1) + 2*cos(3*y+2) + ...
            3*cos(4*y+3) + 4*cos(5*y+4) + 5*cos(4*y+5)) .* ...
            (cos(2*z+1) + 2*cos(3*z+2)+ 3*cos(4*z+3) + 4*cos(5*z+4) + ...
            5*cos(4*z+5));
        f = chebfun3(fa);
        % Try e.g. plot(f) or scan(f)

    case 'rosenbrock'
        % Rosenbrock Function
        % http://www.sfu.ca/~ssurjano/rosen.html
        fa = @(x,y,z) 100*(y-x.^2).^2 + (x-1).^2 + 100*(z-y.^2).^2 + (y-1).^2;
        %f = chebfun3(fa);
        f = chebfun3(fa, [-2 2 -2 2 -2 2]);
        % Try e.g., scan(f)

    case 'aurentz'
        % Nice function suggested by Jared Aurentz. This is a 3D radial 
        % bump function, i.e., a 2D bump which travels around a circle of 
        % radius r about the (x,y)=(0,0).        
        r = 0.6; 
        fa = @(x,y,t) exp(-((x-r*cos(pi*t)).^2 + (y-r*sin(pi*t)).^2)./0.1);
        f = chebfun3(fa);
        % Try e.g., scan(f), or surf(f)  
        
    case 'wave1'
        % A standing wave function
        fa = @(x,y,z) 3*sin(2*x) .* sin(3*y) .* sin(z);
        f = chebfun3(fa);
        
    case 'wave2'
        % A travelling wave function
        fa = @(x,y,t) 3*sin(t + 2*x - 3*y);
        f = chebfun3(fa);
        
    case 'slevinsky'
        % Function suggested by Mikael Slevinsky related to the heat and
        % wave equations
        fa = @(x,y,z) besselj(0,10*sqrt(x.^2+y.^2+z.^2));
        f = chebfun3(fa);
        % Try e.g., surf(f)
        
    case 'levy'
        % http://www.sfu.ca/~ssurjano/levy.html
        fa = @(x,y,z) (sin(pi*(1 + (x - 1)/4))).^2 + ((z - 1)/4).^2 .* ...
            (1+(sin(2*pi*(z - 1)/4)).^2) + (((x - 1)/4)).^2 .* ...
            (1+10*(sin(pi*(1 + (x - 1)/4)+1)).^2) + ...
            (((y - 1)/4)).^2 .* (1+10*(sin(pi*(1 + (y - 1)/4)+1)).^2);
        f = chebfun3(fa, [-10 10 -10 10 -10 10]);
        % Try e.g. surf(f)

       case 'rastrigin'
           % Rastrigin Function (http://www.sfu.ca/~ssurjano/rastr.html)
           fa = @(x,y,z) 30 + (x.^2-10*cos(2*pi*x)) + ...
               (y.^2-10*cos(2*pi*y)) + (z.^2-10*cos(2*pi*z));
           f = chebfun3(fa);

    case 'clebsch'
        % Clebsch Diagonal Cubic from algebraic geometry:
        % http://mathworld.wolfram.com/ClebschDiagonalCubic.html
        % http://blogs.ams.org/visualinsight/2016/03/01/clebsch-surface/
        fa = @(x,y,z) 81*(x.^3 + y.^3 + z.^3) - 189*(x.^2.*y + x.^2.*z +...
            y.^2.*x + y.^2.*z + z.^2.*x + z.^2.*y) + 54*(x.*y.*z) + ...
            126*(x.*y + x.*z + y.*z) - 9*(x.^2 + y.^2 + z.^2) - ...
            9*(x + y + z) + 1;
        f = chebfun3(fa);
        % Try e.g. isosurface(f,0)
        
    case 'cassini'
        % Cassini Surface:
        % https://archive.lib.msu.edu/crcmath/math/math/c/c086.htm
        a = 0.4; 
        fa = @(x,y,z) ((x-a).^2 + y.^2) .* ((x+a).^2 + y.^2) - z.^4;
        f = chebfun3(fa);
        % Try e.g. isosurface(f,0) and isosurface(f,0.5)
        
    case 'barth10'
        % Barth Decic
        % https://archive.lib.msu.edu/crcmath/math/math/b/b054.htm
        a = 1; 
        p = (sqrt(5) + 1) / 2; % the golden section
        fa = @(x,y,z) 8*(x.^2 - p^4*y.^2) .* (y.^2 - p^4*z.^2) .* ...
            (z.^2 - p^4.*x.^2) .* (x.^4 + y.^4 + z.^4 - 2*x.^2.*y.^2 -...
            2*x.^2.*z.^2 - 2*y.^2.*z.^2) + a*(3 + 5*p) .* ((x.^2 + y.^2 ...
            + z.^2 - a)).^2 .* ((x.^2 + y.^2 + z.^2 - (2-p)*a)).^2;
        f = chebfun3(fa);
        % Try e.g. isosurface(f,0)           
           
    case '3dxy'
        % http://csi.chemie.tu-darmstadt.de/ak/immel/script/redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/immel/tutorials/orbitals/index.html
        const = (1/sqrt(4*pi))*sqrt(60/4)*(1*(9*sqrt(30)));
        z0 = 31;
        %z0 = 2;
        fa = @(x,y,z) const * 4./9 * x.*y .* exp(-sqrt(x.^2 + y.^2 + (z-z0).^2)./3);
        %f = chebfun3(fa);
        f = chebfun3(fa, [-30 30 -30 30 -30 30]);
        % Try e.g. isosurface(f) or isosurface(f, [+0.0005 -0.0005])
        
    case '4px'
        % http://csi.chemie.tu-darmstadt.de/ak/immel/script/redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/immel/tutorials/orbitals/index.html
        const = sqrt(3/(4*pi))*(1/(64*sqrt(15)));
        z0 = 31;
        fa = @(x,y,z) const * x .* (20 - 5 * sqrt(x.^2 + y.^2 + (z-z0).^2) + ...
            ((x.^2 + y.^2 + (z-z0).^2))/4) .* ...
                    exp(-sqrt(x.^2 + y.^2 + (z-z0).^2)./4);
        f = chebfun3(fa, [-30 30 -30 30 -30 30]);
        % Try e.g. isosurface(f, [0.0005 -0.0005])
    
    case '32m1'
        % psi_{3,2,-1}, i.e., n = 3, l = 2, m = -1:
        % https://en.wikiversity.org/wiki/Quantum_mechanics/The_hydrogen_atom#Eigenvalues_and_eigenfunctions_of_the_hydrogen_atom
        %psi32_1 = @(r, theta, phi) const * r.^2  .* exp(-r/(3*a0)) .* 
        %                           sin(theta) .* cos(theta) .* exp(-i*phi)
        a0 = 0.0529; % The Bohr radius in nanometer
        %z0 = 31;
        z0 = 2;
        const = 2*sqrt(2)/(81*sqrt(15*a0^7))*sqrt(15/(8*pi));
        fa = @(x, y, z) const * z .* (x + i*y) .* ...
            exp(-sqrt(x.^2 + y.^2 + (z-z0).^2)/(3*a0));
        f = chebfun3(fa);
        %f = chebfun3(fa, [-30 30 -30 30 -30 30]);
        % Try e.g. isosurface(f, [-0.005 0.005])
        
    case '5gz3y'
        % Cartesian Wave Function - 5gz3y
        % http://csi.chemie.tu-darmstadt.de/ak/immel/script/redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/immel/tutorials/orbitals/index.html
        n = 5; 
        z0 = 31;
        const = 1/(900*sqrt(70));
        fa = @(x,y,z) const * y.*z .* (4*(z-z0).^2 - 3*x.^2 - 3*y.^2)...
            ./sqrt(x.^2 + y.^2 + (z-z0).^2).^4 .* ...
            (2./n .* sqrt(x.^2 + y.^2 + (z-z0).^2)).^4 .* ...
            exp(-2/n*sqrt(x.^2 + y.^2 + (z-z0).^2)./2);
        f = chebfun3(fa, [-30 30 -30 30 -30 30]);
        % Try e.g., plot(f) or isosurface(f, [-0.005 0.005])
        
    % Throw an error if the input is unknown.
    otherwise
        %error('CHEB:GALLERY3:unknown:unknownFunction', ...
        error('CHEBFUN:GALLERY3:unknown:unknownFunction', ...
            'Unknown function.')
end

% Only return something if there is an output argument.
if ( nargout > 0 )
    varargout = {f, fa};
else
    % Otherwise, plot the function.
    plot(f)
end

end

%% TODO: Add the following functions:

% 3D Franke
% animation(chebfun3s(@(x,y,z) 3/4*exp(-( (9*x-2).^2 + (9*y-2).^2 + (9*z-2).^2)/4) + 3/4*exp( -((9*x+1).^2)/49 - (9*y+1)/10 - (9*z+1)/10)...
%     +1/2*exp(-( (9*x-7).^2 + (9*y-3).^2 + (9*z-5).^2)/4) - 1/5*exp( -(9*x-4).^2 - (9*y-7).^2 - (9*z-5).^2)))

%%close all; animation(chebfun3s(@(x,y,z) x.^2 .* y.^2 .* z.^5),3) % (https://www.youtube.com/watch?v=w4ixBow28nc)


% close all; animation(chebfun3s(@(x,y,z) exp(-10*(x.^2 + y.^2 + z.^2)) + exp(-2*((x-0.4).^2 + (y-0.3).^2 + (z-0.7).^2)) + exp(-20*((x+0.4).^2 + (y+0.8).^2 + (z-0.85).^2) )))

%http://www.sfu.ca/~ssurjano/easom.html
%ff = @(x,y,z) -cos(x).*cos(y).*cos(z).* exp(-(x-0.5).^2 -(y-0.5).^2 -(z-0.5).^2);

% http://www.sfu.ca/~ssurjano/michal.html
% m = 10; ff = @(x,y,z) -(sin(x).*(sin(x.^2)).^(2*m)) -(sin(y).*(sin(y.^2)).^(2*m)) -(sin(z).*(sin(z.^2)).^(2*m));

%http://www.sfu.ca/~ssurjano/stybtang.html

%http://www.sfu.ca/~ssurjano/detpep10curv.html
%ff = @(x,y,z) 4 * (x - 2 + 8*y - 8*y.^2).^2 + (3 - 4*y).^2 + 16 * sqrt(z+1) .* (2*z-1).^2;

%close all; animation(chebfun3s(@(x,y,z) sin(4*(x.^2 + y.^2 + z.^2))))

%%
%http://www.sfu.ca/~ssurjano/hart3.html
% function y = hart3(x)
% % 
% % Hartmann function 
% % Matlab Code by A. Hedar (Sep. 29, 2005).
% % The number of variables n = 3.
% % 
% a(:,2)=10.0*ones(4,1);
% for j=1:2;
%    a(2*j-1,1)=3.0; a(2*j,1)=0.1; 
%    a(2*j-1,3)=30.0; a(2*j,3)=35.0; 
% end
% c(1)=1.0;c(2)=1.2;c(3)=3.0;c(4)=3.2;
% p(1,1)=0.36890;p(1,2)=0.11700;p(1,3)=0.26730;
% p(2,1)=0.46990;p(2,2)=0.43870;p(2,3)=0.74700;
% p(3,1)=0.10910;p(3,2)=0.87320;p(3,3)=0.55470;
% p(4,1)=0.03815;p(4,2)=0.57430;p(4,3)=0.88280;
% s = 0;
% for i=1:4;
%    sm=0;
%    for j=1:3;
%       sm=sm+a(i,j)*(x(j)-p(i,j))^2;
%    end
%    s=s+c(i)*exp(-sm);
% end
% y = -s;
%%

%% [TODO]: 
% 1) Spherical harmonics: Choose some which have all 3 of theta, phi and r
% and not only theta and phi.
% https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
% Some plots are here:
% https://en.wikipedia.org/wiki/Spherical_harmonics
% As Jared suggests a reasonable way to plot these complex-valued functions 
% Y_L^M is to plot abs((Y_L^M).^2) to get rid of complex values and also 
% plot a physically meaningful concept. 
% In this context we have x = r sin(theta) cos(phi), 
%                         y = r sin(theta) sin(phi), and 
%                         z = r cos(theta).
% A good work to do is to first write an m-file which creates a function 
% handle in terms of X, Y and Z for each Y_L^M (like output of gallery)
% and then call chebfun3 on them.

% 2) Hydrogen atom eigenfunctions as suggested by Mikael
% https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Three-dimensional_examples#Hydrogen_atom

% Each wave function can be factorized into the radial part R(r) and an 
% angular part Y(theta, phi) which is a spherical harmonic, i.e.,
%                                 psi(r, theta, phi) = R(r) * Y(theta, phi)
%
% The eigenvalue problem associated with the hydrogen atom is the following
%               \hat{H} psi(r, theta, phi) = E psi(r, theta, phi)
% where \hat{H} is the Hamiltonian ...
% E := E_n = -0.5 \alpha^2 \mu c^2 / n^2
% \alpha = fine structure constant = 1/137

% a0 = 0.529; % The Bohr radius in angstroms
%a0 = 0.529e-10; % The Bohr radius in meters

% Each wave function has three quantum numbers:
% n = principal quantum number (quantum number for the radial part r)
% l = azimuthal quantum number
% m = magnetic quantum number

% nucleeus (i.e., proton in the case of Hydrogen atom) is fixed at the 
% origin (x0, y0, z0), because it is much more heavier.
% r = sqrt(x.^2 + y.^2 + z.^2); distance of the electron from the nucleus
% (i.e, proton in case of Hydrogen atom)
% theta = elevation
% phi = azimuthal angle

% x = r * sin(theta) * cos(phi)
% y = r * sin(theta) * sin(phi)
% z = r * cos(phi)
% tan(theta) = sqrt(x^2+y^2) / z
% tan(phi) = y / x
%%

