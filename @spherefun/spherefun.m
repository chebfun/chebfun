classdef spherefun < separableApprox
    
    % TODO: Improve documentation of input options.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = spherefun(varargin)
            % The main spherefun constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                f.domain = [-pi pi 0 pi];
                return
            end
            
            % Type of construction
            constructorType = 1; % Default using rank BMC preserving rank 1 updates.
            
            % Remove this code when we are done testing the constructor.
            if numel(varargin) > 1
                if ischar(varargin{2})
                    if strcmpi(varargin{2},'2by2')
                        constructorType = 2;
                        varargin(2) = [];
                    elseif strcmpi(varargin{2},'PlusMinus')
                        constructorType = 3;
                        varargin(2) = [];
                    end
                end
            end
            
            % Call the constructor, all the work is done here:
            switch constructorType
                case 1
                    f = constructor(f, varargin{:});
                case 2
                    f = constructor2by2(f, varargin{:});
                case 3
                    f = constructorPlusMinus(f, varargin{:});
                otherwise
                    f = constructor(f, varargin{:});
            end
                    
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
    methods ( Access = public, Static = false, Hidden = true )
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % The main bulk of the SPHEREFUN constructor:
        g = constructor(g, op, dom, varargin);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Poisson solver: 
        u = Poisson( f, const, n);
        
        % Converts a function in spherical coordinates to one in Cartesian
        % coordinates on the sphere.
        fdf = sphf2cartf(f, lam, th, coord);
        
        % Degree l Order m spherical harmonic.
        Y = sphharm(l,m,coord)          
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by SPHEREFUN class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % DOMAIN: default for the is [-pi,pi] x [0,pi].
        % which corresponds to using colatitude for the elevation angle
        % (second input argument).  Doubled-up sphere will have a domain of
        % [-pi,pi] x [-pi,pi].        
        idxPlus
        idxMinus
        pivotIndices
        nonZeroPoles = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private constant properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Constant )
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