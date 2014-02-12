% CHEBFUN2V Class constructor for CHEBFUN2V objects
% 
% CHEBFUN2V(F,G) constructs a CHEBFUN2V with two components from the function handles F
% and G.  F and G can also be CHEBFUN2 objects or any other object that the
% CHEBFUN2 constructor accepts.  Each component is represented as a CHEBFUN2. 
%
% CHEBFUN2V(F,G,H) constructs a CHEBFUN2V with three components from the
% function handles F, G, and H.  F, G, and H can also be CHEBFUN2 objects 
% or any other object that the CHEBFUN2 constructor accepts. 
%
% CHEBFUN2V(F,G,[A B C D]) constructs a CHEBFUN2V object from F and G 
% on the domain [A B] x [C D].
%
% CHEBFUN2V(F,G,H,[A B C D]) constructs a CHEBFUN2V object from F, G, and 
% H on the domain [A B] x [C D].
% 
% See also CHEBFUN2. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

classdef chebfun2v
    
    properties ( Access = public )
        components   % Array of CHEBFUN2 objects. 
        nComponents  % Number of components
        isTransposed % transposed?
    end
    
    methods
        
        function f = chebfun2v(varargin)
            % The main CHEBFUN2V constructor!
            
            % Return an empty CHEBFUN2V:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % constructor.m does all the work: 
            f = constructor(f, varargin{:});
            
        end
    end 
    
end