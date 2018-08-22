classdef ballfun
%BALLFUN   BALLFUN class for representing scalar-valued functions on the unit ball.
%
%  Class for approximating smooth scalar-valued functions defined on the unit ball.
%
% BALLFUN(F) constructs a BALLFUN object representing the scalar-valued function F
% on the unit ball. F may be a function handle or scalar. If F is a function handle,
% then it should be vectorized.
%
% See also SPHEREFUN, DISKFUN, BALLFUNV.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASS PROPERTIES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        coeffs   % The Chebyshev-Fourier-Fourier coefficients of a BALLFUN
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = ballfun(X,varargin)
            % The constructor.

            if (nargin == 1) && isa(X, 'double')
                % X = tensor of Chebyshev-Fourier-Fourier coefficients. in
                % r, lam, th, where lam is the azimuthal variable and th the polar angle

                f.coeffs = X;

            elseif (nargin == 2) && isa(X, 'function_handle') && (nargin(X)==3) && isa(varargin{1}, 'double')
                % X = function_handle of r,lam,th, varargin{1} = [m,n,p]
                g = ballfun.fun2ballfun(X,varargin{1});

                f.coeffs = g.coeffs;

            elseif (nargin == 3) && isa(X, 'function_handle') && (nargin(X)==3) && strcmp(varargin{1}, 'cart') && isa(varargin{2}, 'double')
                % X = function_handle of x,y,z, varargin{1} = 'cart', varargin{2} = [m,n,p]

                x = @(r,lam,th)r.*sin(th).*cos(lam);
                y = @(r,lam,th)r.*sin(th).*sin(lam);
                z = @(r,lam,th)r.*cos(th);
                g = ballfun.fun2ballfun(@(r,lam,th)X(x(r,lam,th),y(r,lam,th),z(r,lam,th)),varargin{2});
                f.coeffs = g.coeffs;

            elseif (nargin == 2) && isa(X, 'double') && (varargin{1} == "vals")
                % X = tensor of values, varargin{1} = 'vals' ; X(r,lam,th),
                % -1 <= r <= 1, -pi <= lam <= pi, -pi <= th <= pi,

                f.coeffs = ballfun.vals2coeffs(X);
            else

                error('BALLFUN:BALLFUN:unknown', ...
                ['Undefined function ''ballfun'' for input arguments of type ' ...
                '%s.'], class(X));

            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATIC METHODS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Convert to Chebyshev--Fourier--Fourier values
        VALS = coeffs2vals(CFS);
        
        % Convert to Chebyshev--Fourier--Fourier values
        CFS = vals2coeffs(VALS);
        
        % Convert a chebfun function to a ballfun function
        f = chebfun2ballfun(g,S);
        
        % Convert a spherefun function to a ballfun function
        f = spherefun2ballfun(g,S);
        
        % Compute the whole set of normalized associated Legendre
        % polynomials
        A = normalized_legendre(n);
        
        f = solHarm(l,m);
        
    end
    
    methods ( Access = private, Static = true )
        
        % Convert a handle function to a ballfun function
        f = fun2ballfun(g,S);
        
    end
    
    methods ( Access = private, Static = false )
        
    end
    
end
