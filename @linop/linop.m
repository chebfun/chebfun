classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) linop < chebmatrix
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    
    properties
        constraint = linopConstraint()
        continuity = linopConstraint()
    end
    
    methods
        function L = linop(M)
            %LINOP   Linop constructor.
            %   Linops are typically onstructed from a LINBLOCK or a CHEBMATRIX.
                
            % TODO: Better argument checking
            L = L@chebmatrix(M);
            
         end
                
    end
    
    % These are provided as more convenient names than the linBlock equivalents.
    methods (Static)
        function D = diff(varargin)
            D = linop( linBlock.diff(varargin{:}) );
        end
        
        function C = cumsum(varargin)
            C = linop( linBlock.cumsum(varargin{:}) );
        end
        
        function I = eye(varargin)
            I = linop( linBlock.eye(varargin{:}) );
        end
        
        function Z = zeros(varargin)
            Z = linop( linBlock.zeros(varargin{:}) );
        end
        
        function U = mult(varargin)
            U = linop( linBlock.mult(varargin{:}) );
        end
        
        function Z = zero(varargin)
            Z = linop( linBlock.zero(varargin{:}) );
        end
        
        function S = sum(varargin)
            S = linop( linBlock.sum(varargin{:}) );
        end
        
        function E = feval(varargin)
            E = linop( linBlock.feval(varargin{:}) );
        end
        
        function E = eval(varargin)
            E = @(x) linop.feval(x, varargin{:});
        end
        
        function F = inner(varargin)
            F = linop( linBlock.inner(varargin{:}) );
        end
        
        function F = dot(varargin)   % synonym for inner()
            F = linop( linBlock.dot(varargin{:}) );
        end
        
        function F = fred(varargin)   % synonym for inner()
            F = linop( linBlock.fred(varargin{:}) );
        end

        function V = volt(varargin)   % synonym for inner()
            V = linop( linBlock.volt(varargin{:}) );
        end
        
    end
    
end