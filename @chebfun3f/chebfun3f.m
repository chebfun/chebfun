classdef chebfun3f

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASS PROPERTIES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties   
    % COLS: Mode-1 fibers, i.e., columns which are functions of x used in 
    % Tucker representation.
    cols
    
    % ROWS: Mode-2 fibers, i.e. rows which are functions of y used in 
    % Tucker representation.
    rows
    
    % TUBES: Mode-3 fibers, i.e. tubes which are functions of z used in 
    % Tucker representation.
    tubes
    
    % CORE: discrete core tensor in Tucker representation
    core
    
    % DOMAIN: box of CHEBFUN3, default is [-1, 1] x [-1, 1] x [-1, 1].
    domain
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASS CONSTRUCTOR:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    function cf3f = chebfun3f(varargin)
        cf3f = constructor(cf3f, varargin{:});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASS METHODS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access = public)
    varargout = rank(cf3F);
    varargout = degree(cf3F);
    varargout = feval(cf3F,x,y,z);
end

end