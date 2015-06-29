classdef (InferiorClasses = {?chebfun}) domain < double
%DOMAIN   Utility class for CHEBFUN. Mostly for backward compatibility.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOMAIN Class Description:
%
% DOMAIN inherits from a standard Matlab DOUBLE. A domain object only
% contains vector for the endpoints and breakpoints of the interval it
% represents. This class is lightly documented, since it is mostly intended
% for backward compatibility.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function obj = domain(varargin)
            %Constructor for the DOMIAN class.
            
            % Return an empty DOMAIN on null input:
            if ( nargin == 0 )
                data = [];
            else
                data = horzcat(varargin{:});
            end
            
            % Create the domain:
            obj = obj@double(data);
            
        end
        
    end     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function out = validate(dom)
            out = true;
            if ( ~isnumeric(dom) || any(isnan(dom)) )
                out = false;
            elseif ( isempty(dom) )
                return
            elseif ( size(dom, 1) > 1 || size(dom, 2) < 2 )
                out = false;
            elseif ( any(diff(dom) <= 0) )
                out = false;
            end
            if ( out == false && nargout == 0 )
                error('CHEBFUN:DOMAIN:domain:invalid', ...
                    'Ends must be a 1xM vector of ordered doubles.');
            end
        end
        
        function display(dom)
            disp(double(dom))
        end
        
        function out = subsref(d, s)
            out = subsref@double(d, s);
            out = double(out);
        end
        
        function varargout = sprintf(varargin)
            % This is required as built-in subsref does not know what to do with
            % a DOMAIN object.
            
            varargin = domain2double(varargin{:});
            
            % Call built-in SPRINTF:
            varargout{1:nargout} = sprintf(varargin{:});
        end
                   
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Merge two domains:
        newDom = merge(varargin)
        
        function varargin = toDouble(varargin)
            % Cast DOMAIN to DOUBLE:
            for k = 1:nargin
                if ( isa(varargin{k}, 'domain') )
                    varargin{k} = double(varargin{k});
                end
            end
        end
        
    end
    
end
