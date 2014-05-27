classdef(InferiorClasses = {?chebfun2}) ADchebfun2
    %ADCHEBFUN2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        chebfun2;   % ADchebfun2 has a chebfun2
        der;
    end
    
    methods
        % Main constructor. Convert a chebfun2 to ADchebfun2
        function g = ADchebfun2 ( varargin )
            if( nargin == 0 )
                % return an empty chebfun2 object. 
            elseif isa(varargin{1},'chebfun2')
                g.chebfun2 = varargin{:};  % Assign to the chebfun2 field of g.
                g.der = chebfun2der(1,g.chebfun2.domain);
            else
                cheb2temp = chebfun2(varargin{:});
                g = ADchebfun2(cheb2temp);
            end
        end
        
        function varargout = plot(f)
            if nargout
                varargout = plot(f.chebfun2);
            else
                plot(f.chebfun2);
            end
        end
    end
    
end

