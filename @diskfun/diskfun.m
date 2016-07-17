classdef diskfun < separableApprox
% DISKFUN class for representing functions on the unit disk.
% 
%   Class for approximating functions defined on the unit disk. The 
%   functions should be smooth.
%
% DISKFUN(F) constructs a DISKFUN object representing the function F on
% the unit. F can have the following form:
%    1. A function handle in Cartesian coordinates (x,y), e.g., @(x,y) x.*y + cos(x).
%    2. A function handle in polar coordinates (theta,r), where
%       theta is the angular variable and satisfies -pi <= theta <= pi
%       and r is the radial variable and satisfies 0 <= r < 1,
%       e.g., @(theta,r) cos(r.*sin(theta))
%    3. A matrix of numbers. 
% If F is a function handle then it should allow for vectorized 
% evaluations.
%
% If F is a matrix, F = (f_ij), the numbers fij are used as function values
% at tensor Fourier-Chebyshev points in the polar coordinate 
% system (theta, r) in [-pi,pi]x[0,1].
%
% The DISKFUN software system is based on: 
%
% A. Townsend, H. Wilber, and G. Wright, Computing with function on
% spherical and polar geometries II: The disk, submitted, 2016. 
%
% See also CHEBFUN2, SPHEREFUN, DISKFUNV.

    
    % TODO: Improve documentation of input options.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = diskfun(varargin)
            % The main diskfun constructor!
            
            % Return an empty DISKFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                f.domain = [-pi pi 0 1];
                return
            end
            
            % Type of construction
            constructorType = 1; % Default using rank BMC preserving rank 1 updates.
            
            % Remove this code when we are done testing the constructor.
            if numel(varargin) > 1
                if ischar(varargin{2})
                    if strcmpi(varargin{2},'old')
                        constructorType = 2;
                        varargin{2} = [];
                    end
                end
            end
            
            % Call the constructor, all the work is done here:
            switch constructorType
                case 1
                    f = constructor(f, varargin{:});
                case 2
                    f = constructoroldII(f, varargin{:});
                otherwise
                    f = constructor(f, varargin{:});
            end
                    
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function f = set.coords(f, propName)
            if (strcmp(propName, 'polar') || strcmp(propName, 'cart'))
                f.coords = propName;
            else  %error if unacceptable setting is provided
    error('CHEBFUN:DISKFUN:setcoords:propName', ...
            'Coordinate setting must be either ''polar'' or ''cart''')
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        
    end
      
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
     
        % The main bulk of the DISKFUN constructor:
        g = constructor(g, op, dom, varargin);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Fast Poisson solver: 
        u = poisson( f, bc, m, n );
        
        % Converts a function in polar coordinates to one in Cartesian
        % coordinates on the disk
        fdf = pol2cartf(f, r, th);
        
        % Disk harmonics
        Y = diskharm(l,m,type) ;   
        
        %coordsettings: determines whether a function has input/wants output in
         %polar or cartesian coords. 
        
        iscart = coordsetting(varargin);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by DISKFUN class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
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
        coords = 'cart'; %default: results returned in Cartesian coordinates
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private constant properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Constant )
        %alpha = 50;  % Growth factor control.
    end
    
   
end