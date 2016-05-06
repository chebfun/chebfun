classdef spherefun < separableApprox
%SPHEREFUN class for representing functions on the unit sphere.
% 
%   Class for approximating functions defined on the unit sphere. The 
%   functions should be smooth.
%
% SPHEREFUN(F) constructs a SPHEREFUN object representing the function F on
% the unit. F can have the following form:
%    1. A function handle in (x,y,z), e.g., @(x,y,z) x.*y.*z + cos(x).
%    2. A function handle in spherical coordinates (lambdda,theta), where
%       lambda is the azimuthal variable and satisfies -pi <= lambda <= pi
%       and theta is the polar angle and satisfies 0 <= theta < pi,
%       e.g., @(lambda,theta) cos(cos(lambda).*sin(theta))
%    3. A matrix of numbers. 
% If F is a function handle then it should allow for vectorized 
% evaluations.
%
% If F is a matrix, F = (f_ij), the numbers fij are used as function values
% at tensor equally-spaced points in the intrinsic spherical coordinate 
% system, i.e., [-pi,pi]x[0,pi].
%
% The SPHEREFUN software system is based on: 
%
% A. Townsend, H. Wilber, and G. Wright, Computing with function on
% spherical and polar geometries I: The sphere, submitted, 2015. 
%
% See also CHEBFUN2, SPHEREFUNV.

% Copyright 2016 by The University of Oxford and The CHEBFUN Developers.
% See http://www.chebfun.org/ for CHEBFUN information.

% TODO: Improve documentation of input options. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = false)
        
        function f = spherefun(varargin)
            % The main spherefun constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                f.domain = [-pi pi 0 pi];
                return
            end
            
            f = constructor(f, varargin{:});                    
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %f = conj(f);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = false, Hidden = true)
       f = projectOntoBMCI(f);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = false)
        
        % The main bulk of the SPHEREFUN constructor:
        g = constructor(g, op, dom, varargin);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = true)
  
        % Poisson solver: 
        u = poisson(f, const, m, n);
        
        % Helmholtz solver: 
        u = helmholtz(f, K, m, n);
        
        % Convert matrix of coefficients to a spherefun: 
        f = coeffs2spherefun(CFS); 
        
        % Convert a function in spherical coordinates to one in Cartesian
        % coordinates on the sphere.
        fdf = sphf2cartf(f, lam, th, coord);
        
        % Degree l Order m spherical harmonic.
        Y = sphharm(l, m, coord);
        
        % Plot the outline of the landmasses of earth
        h = plotEarth(linespec);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by SPHEREFUN class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private, Static = true)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        % DOMAIN: default is [-pi,pi] x [0,pi] which corresponds to using 
        % colatitude for the elevation angle (second input argument). 
        % Doubled-up sphere will have a domain of [-pi,pi] x [-pi,pi].
        idxPlus
        idxMinus
        nonZeroPoles = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private constant properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Constant)
%         alpha = 2;  % Growth factor control.
    end
    
    %
    %     methods
    %         function g = spherefun( varargin )
    %             if( nargin == 0 )
    %
    %             else
    %                 g = constructor(g , varargin{:} );  % pass to constructor.
    %             end
    %         end
    %     end
    
end