function varargout = gallery3(name)
%CHEB.GALLERY3   CHEBFUN3 example functions.
%
%   F = CHEB.GALLERY3(NAME) returns a chebfun3 corresponding to NAME.
%   The domain is often but not always [-1,1]^3.
%   See the listing below for available names.
%
%   [F, FA] = CHEB.GALLERY3(NAME) also returns the anonymous function FA 
%   used to define the chebfun3.
%
%   CHEB.GALLERY3 with no output argument creates a plot of the selected
%   function. 
%
%   Below, the functions are listed with suggestions for good
%   viewing in the form e.g. [plot(f)].  For more details about
%   the functions, see the code and associated comments.
%
%   aurentz      A 3D radial bump function [scan(f)]
%   barth10      Barth decic from algebraic geometry
%   bessel       Bessel function of radius r
%   cassini      Cassini Surface from algebraic geometry
%   clebsch      Clebsch diagonal cubic from alg. geom. [isosurface(f,0)]
%   genz1        Oscillatory integration test fun. [plot(f)]
%   genz2        Cubature test function with corner peak [surf(f)]
%   genz3        Cubature test function with Gaussian peak [surf(f)]
%   genz4        Cubature test function with product peak [surf(f)]
%   griewank     Function with many local minima.  The global
%                min. is 0 at (0,0,0). [plot(f)]
%   levy         Function with global min. 0 at (1,1,1) [surf(f)]
%   rastrigin    Function with many local minima.  The global
%                min. is 0 at (0,0,0).  [isosurface(f)]
%   rosenbrock   Function with global min 0 at (1,1,1).
%                [contour(f(:,:,1),0:10:100)]
%   runge        3D Runge function [slice(f)]
%   shubert      Function with many local minima.  [scan(f)]
%   wagon        3D generalization of "SIAM 100 Digit Challenge" function
%                [contourf(f(-.16,:,:))]
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
        'genz4', 'griewank', 'rastrigin', 'bessel'};
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
        % This function was proposed by Stan Wagon.
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
        
    case 'bessel'
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
    title(name)
end

end

