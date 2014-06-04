classdef (InferiorClasses = {?chebfun2}) adchebfun2
%ADCHEBFUN2   A class for supporting automatic differentiation in Chebfun2.
%
%   The ADCHEBFUN2 class allows chebop2 to compute variable coefficients of 
%   partial differential operators.
% 
%   This class is not intended to be called directly by the end user.
%
%   V = ADCHEBFUN2(U), where U is a CHEBFUN2, returns the ADCHEBFUN2
%   object V, which has a derivatives seeded as the identity
%   operator on the domain of U.
%
%   V = ADCHEBFUN2(...), where the input is the same as would be
%   passed to the CHEBFUN2 constructor, constructs a CHEBFUN2 from the
%   input to the method. It then returns the ADCHEBFUN object V with
%   the function part consisting of the CHEBFUN constructed, and the
%   derivative part as the identity operator on the domain.
%
%   THIS CLASS HAS LIMITED FUNCTIONALITY AND IS DESIGNED SPECIFICALLY FOR 
%   CHEBOP2 REQUIREMENTS ONLY. 
%
% See also ADCHEBFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    properties
        chebfun2   % An adchebfun2 has a chebfun2
        der        % Derivative information stored here
    end
    
    methods
        % Main constructor. Convert a chebfun2 to ADchebfun2
        function g = adchebfun2 ( varargin )
            if( nargin == 0 )
                % return an empty chebfun2 object. 
            elseif isa(varargin{1},'chebfun2')
                g.chebfun2 = varargin{:};  % Assign to the chebfun2 field of g.
                g.der = chebfun2der(1,g.chebfun2.domain);
            else
                cheb2temp = chebfun2(varargin{:});
                g = adchebfun2(cheb2temp);
            end
        end
        
        function varargout = plot(f)
            % Plot a chebfun2 
            if nargout
                varargout = plot(f.chebfun2);
            else
                plot(f.chebfun2);
            end
        end
    end
    
end

