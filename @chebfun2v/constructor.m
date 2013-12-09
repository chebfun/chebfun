function G = constructor(G, op, varargin)

if ( nargin == 3 )
    domain = varargin{1}; 
else
    domain = [-1 1 -1 1];
end
G.components = []; 
for j = 1:numel( op )
    f = chebfun2(op{ j }, domain);  
    G.components = [G.components f]; 
end
G.nComponents = numel( op ); 

end