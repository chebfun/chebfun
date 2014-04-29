function G = constructor(G, op, varargin)
% CTOR  CHEBFUN2V constructor
% This function calls the CHEBFUN2 constructor once for each non-zero
% component because a CHEBFUN2V is just vector of CHEBFUN2 objects.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if isa(op,'chebfun2v') % argument is a CHEBFUN2V so there is nothing to do.
    G = op;
    return
end

domain = [-1 1 -1 1];

if ( nargin >= 3 )
    if ( isa(varargin{1}, 'double') )
        % Grab domain. 
        domain = varargin{1};
    elseif (isa( op, 'function_handle' ) && (isa( varargin{1}, 'function_handle' ) ) )
        % The constructor is given the components as function handles. 
        len = numel(varargin); 
        j = 0; 
        while ( ( j < len) && isa( varargin{j+1}, 'function_handle' ) )
             j = j + 1;
        end
        op = {op, varargin{1:j}};  
        if (j < len)
            domain = varargin{j+1}; 
        end
        if ( j + 1 < len )
            varargin = varargin{j+2:len}; 
        else
            varargin = []; 
        end
    elseif ( isa( op, 'chebfun2' ) && isa( varargin{1}, 'chebfun2' ) ) 
        % The constructor is given the components as CHEBFUN2 objects. 
        domain = op.domain; 
        len = numel(varargin); 
        j = 0; 
        while ( ( j < len) && isa( varargin{j+1}, 'chebfun2' ) )
            if ( ~domainCheck(op, varargin{j+1}) )
                error('CHEBFUN2V:CONSTRUCTOR','Domain not the same.');
            end
            j = j + 1;
        end 
        op = {op, varargin{1:j}};  
        if (j < len)
            varargin = varargin{j+1:len}; 
        else
            varargin = []; 
        end 
    else
        error('CHEBFUN2V:CONSTRUCTOR','Unrecognised syntax.');
    end
end

% Pull out domain
for j = 1:numel( op ) 
    if ( isa( op{j}, 'chebfun2' ) )
        domain = op{j}.domain; 
    end
end

% Loop through the cell array of components. 
G.components = [];
for j = 1:numel( op )
    f = chebfun2(op{ j }, domain, varargin );
    G.components{j} = f;
end
G.nComponents = numel( op );
G.isTransposed = 0; 

if ( numel(op) > 3) 
    error('CHEBFUN2:CONSTRUCTOR:ARRAYVALUED','More than three components is not supported.')
end
end
