classdef ballfunv
%BALLFUNV   BALLFUNV class for representing vector-valued functions on the unit ball.
%
%  Class for approximating smooth vector-valued functions defined on the unit ball.
%
% BALLFUNV(F, G, H) constructs a BALLFUNV object representing the vector-valued
% function [F;G;H] on the unit ball. F, G, and H may be BALLFUN objects, function
% handles or scalars. If they are function handles, then they those function handles
% should be vectorized.
%
% See also BALLFUN.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASS PROPERTIES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        comp   % The three components Vx, Vy, Vz of a BALLFUNV
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        % Take the 3 components
        function F = ballfunv(varargin)
            
            % Return an empty BALLFUNV:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
            
            % If the argument is a BALLFUNV, nothing to do:
            if ( isa(varargin{1}, 'ballfunv') )
                F = varargin{1};
                return
            end
            
            % BALLFUNV objects are vector-valued so complain if there 
            % are less than 3 components: 
            if ( numel(varargin) < 3 )
                error('BALLFUNV:ballfunv', ...
                    'Less than three components is not supported.')
            end
            
            % BALLFUNV objects cannot contain more than five components. 
            % Complain if we have been given six or more.  
            if ( numel(varargin) > 5 )
                error('BALLFUNV:ballfunv', ...
                    'More than four components is not supported.')
            end
            
            % Create a BALLFUNV from 3 BALLFUN objects
            if ( numel(varargin) == 3 )
                for jj = 1:3
                   if isa(varargin{jj}, 'ballfun') == 0
                       error('BALLFUNV:ballfunv', ...
                        'Input arguments should be ballfun functions')
                   end
                end

                fh{1} = varargin{1};
                fh{2} = varargin{2};
                fh{3} = varargin{3};
            end
           
            % Create a BALLFUNV from 3 function handles and an array
            if ( numel(varargin) == 4 )
                S = varargin{4};
                fh{1} = ballfun(varargin{1},S);
                fh{2} = ballfun(varargin{2},S);
                fh{3} = ballfun(varargin{3},S);
            end
            
            % Create a BALLFUNV from 3 function handles and an array
            if ( numel(varargin) == 5 )
                if ~strcmp(varargin{4}, 'cart')
                    error('BALLFUNV:ballfunv', ...
                        'Input arguments not supported')
                end
                S = varargin{5};
                fh{1} = ballfun(varargin{1},'cart',S);
                fh{2} = ballfun(varargin{2},'cart',S);
                fh{3} = ballfun(varargin{3},'cart',S);
            end
            
            F.comp = fh;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATIC METHODS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        v = PT2ballfunv(P,T);
        v = rcurl(P);
    end
    
    methods ( Access = private, Static = true )
        u = PTequation(f);
    end

    
end
