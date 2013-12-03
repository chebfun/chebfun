classdef linop < chebmatrix
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    
    properties
        constraint = linopConstraint()
        continuity = linopConstraint()
    end
    
    properties (Dependent)
        sizeReduction
    end
    
    methods
        function L = linop(M)
            %LINOP   Linop constructor.
            %   Linops are typically onstructed from a LINBLOCK or a CHEBMATRIX.
                
            % TODO: Better argument checking
            L = L@chebmatrix(M);
            
         end
        
        function s = get.sizeReduction(L)
            s = getDownsampling(L);
        end    
        
    end
    
    % These are provided as more convenient names than the linBlock equivalents.
    methods (Static)
        function D = diff(varargin)
            D = linBlock.diff(varargin{:});
        end
        
        function C = cumsum(varargin)
            C = linBlock.cumsum(varargin{:});
        end
        
        function I = eye(varargin)
            I = linBlock.eye(varargin{:});
        end
        
        function Z = zeros(varargin)
            Z = linBlock.zeros(varargin{:});
        end
        
        function U = mult(varargin)
            U = linBlock.mult(varargin{:});
        end
        
        function Z = zero(varargin)
            Z = linBlock.zero(varargin{:});
        end
        
        function S = sum(varargin)
            S = linBlock.sum(varargin{:});
        end
        
        function E = feval(varargin)
            E = linBlock.feval(varargin{:});
        end
        
        function E = eval(varargin)
            E = linBlock.eval(varargin{:});
        end
        
        function F = inner(varargin)
            F = linBlock.inner(varargin{:});
        end
        
        function F = dot(varargin)   % synonym for inner()
            F = linBlock.dot(varargin{:});
        end
        
    end
    
end