classdef ballfunv
    %BALLFUNV  A library for the DFS method in the solid sphere.
    
    properties
        
        comp   % Components Vx, Vy, Vz of a BALLFUNV
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        % Take the 3 components
        function F = ballfunv(varargin)
            
            % Return an empty BALLFUNVCART:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
            
            % If the argument is a BALLFUNV, nothing to do:
            if ( isa(varargin{1}, 'ballfunv') )
                F = varargin{1};
                return
            end
            
            % BALLFUNVCART objects are vector-valued so complain if there 
            % are less than 3 components: 
            if ( numel(varargin) < 3 )
                error('BALLFUNV:ballfunv', ...
                    'Less than three components is not supported.')
            end
            
            % BALLFUNV objects cannot contain more than three components. 
            % Complain if we have been given four or more.  
            if ( numel(varargin) > 3 )
                error('BALLFUNV:ballfunv', ...
                    'More than three components is not supported.')
            end
            
            for jj = 1:3
               if isa(varargin{jj}, 'ballfun') == 0
                   error('BALLFUNV:ballfunv', ...
                    'Input arguments should be ballfun functions')
               end
            end
            
            fh{1} = varargin{1};
            fh{2} = varargin{2};
            fh{3} = varargin{3};
            F.comp = fh;
        end
    end
    
    methods ( Access = public, Static = true )
        v = PT2ballfunv(P,T);
        v = rcurl(P);
    end
    
    methods ( Access = private, Static = true )
        u = PTequation(f);
    end
    
    methods ( Access = private, Static = false )
        

    end
    
end
