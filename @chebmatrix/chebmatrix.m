classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) chebmatrix
    % No size/compatability checking whatsoever!
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    
    properties
        blocks = {};
        discretizer = @colloc2;        
    end
    
    properties (Dependent)
        domain
        diffOrder
    end
    
    methods
        
        % Constructor.
        function A = chebmatrix(data)
            if isempty(data)
                return
            elseif isa(data,'chebmatrix')
                A.blocks = data.blocks;
            elseif isa(data,'chebfun') || isa(data,'linBlock')
                A.blocks = {data};
            elseif iscell(data)
                A.blocks = data;
            end
            
            % Run this to check domain compatability
            A.domain;
        end
        
        function d = get.domain(L)
            isnum = cellfun(@isnumeric,L.blocks);
            blocks = L.blocks(~isnum);
            d = cellfun(@(x) x.domain,blocks,'uniform',false);
            d = chebfun.mergeDomains(d{:}); 
        end
        
        function d = getDomain(L)
            % DOMAIN(L) returns the domain on which functions are defined for
            % the chebmatrix L.
            d = L.domain;
        end       
        
        function d = get.diffOrder(L)
            d = getDiffOrder(L);
        end
               
        function d = getDiffOrder(A)
            d = zeros(size(A));
            for j = 1:numel(A.blocks);
                if ( isa(A.blocks{j}, 'operatorBlock') )
                    d(j) = A.blocks{j}.diffOrder;
                end
            end
        end
        
        function out = isFunVariable(A,k)
            [rowSize, colSize] = blockSizes(A);
            out = isinf(colSize(1, :));
            if ( nargin > 1 )
                out = out(k);
            end
        end
        
    end
              
end