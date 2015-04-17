classdef spherefun
    
    properties ( GetAccess = 'public' , SetAccess = 'public' )
        Cols            % Cols calculated during GE 
        Rows            % Rows calculated during GE
        BlockDiag       % Pivot matrices used during GE
        PivotLocations  % Locations used during GE
    end
    
    methods 
        function g = spherefun( varargin )
            if( nargin == 0 )
                
            else
                g = constructor(g , varargin{:} );  % pass to constructor. 
            end
        end
    end
    
end