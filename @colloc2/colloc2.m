classdef colloc2 < linopDiscretization
     
    methods
        function disc = colloc2(source)
            if isempty(source)
                return
            end

            % Decide which kind of object we're discretizing:
            if isa(source,'linop')
                disc.linop = source;
            else
                disc.source = source;  % block
            end
            disc.domain = source.domain;
        end
         
        % Specific to blockDiscretization
        D = diff(A, m)
        C = cumsum(A, m)
        F = mult(disc, f)
       
        function I = eye(disc)
            n = disc.dimension;
            I = eye(sum(n));
        end
        
        function Z = zeros(disc)
            n = disc.dimension;
            Z = zeros(sum(n));
        end
        
        function Z = zero(disc)
            n = disc.dimension;
            Z = zeros(1, sum(n));
        end
        
        S = sum(disc)
        
        E = feval(disc, location, direction)
        
        function F = inner(disc, f)
            d = chebmatrix.mergeDomains({disc, f});
            [x, w] = points(disc);
            F = w.*f(x);
        end
        
        
        function C = mtimes(A, B)
            C = A*B;
        end
        
        function C = plus(A, B)
            C = A+B;
        end
        
        function B = uminus(A)
            B = -A;
        end
        
        function fx = toValues(disc,f)
            n = disc.dimension;

            x = points(disc);
            fx = f(x);

            % Evaluate left- and right-sided limits at breaks:
            csn = [0, cumsum(n)];
            dxloc = csn(2:end-1);
            fx(dxloc) = feval(f, x(dxloc), 'left');
            fx(dxloc+1) = feval(f, x(dxloc), 'right');
        end
               
        function L = discretize(disc,dimension)
            A = disc.source;
            if ( nargin > 1 )
                disc.dimension = dimension;
            end
            validate(disc);
            
            if isa(A,'chebmatrix')
                % Evaluate recursively, block by block:
                L = cellfun(@(x) blockDiscretize(disc,x),A.blocks,'uniform',false);
            else
                L = blockDiscretize(disc,A);
            end
        end
                        
        function f = toFunction(disc,values)
            if ( disc.numIntervals > 1 )
                values = mat2cell(values,disc.dimension);
            end
            f = chebfun(values, disc.domain);
        end
               
        % Specific to linopDiscretization        
        function [A,P,B] = matrix(disc)
            if isempty(disc.linop)
                A = discretize(disc);
                return
            end
            L = disc.linop;
            disc.source = L.operator;
            blocks = discretize(disc);
            [rows,P] = disc.reproject(blocks);
            P = blkdiag(P{:});

            B = [];
            if ~isempty(L.constraint)
                disc.source = L.constraint.operator;
                constr = discretize(disc);
                B = [ cell2mat(constr); B ];
            end
            if ~isempty(L.continuity)
                disc.source = L.continuity.operator;
                constr = discretize(disc);
                B = [ cell2mat(constr); B ];
            end
            A = cell2mat(rows);
            A = [ B; A ];

        end
        
        function b = rhs(disc,f)
            if isempty(disc.dimension)
                error('Discretization dimension not given.')
            end
            % NONONO
            if isa(f,'chebfun'), f = chebmatrix({f}); end
            disc.source = f;
            row = discretize(disc);
            if ( ~iscell(row) )
                row = {row};
            end
            row = disc.reproject(row);
            
            b = cell2mat(row);
            L = disc.linop;
            if ~isempty(L.constraint)
                b = [ L.constraint.values; b ];
            end
            if ~isempty(L.continuity)
                b = [ L.continuity.values; b ];
            end
        end
        
    end
    
    methods ( Static )
            
        function [isDone, epsLevel] = testConvergence(v)
            % TODO: (for breakpoints and systems)
            f = chebtech2(v);
            [isDone, epsLevel] = strictCheck(f, 1e-10);
        end
        
    end
    
    methods ( Access = private )
        function [B,P] = reproject(disc,blocks)
            reduce = disc.linop.sizeReduction;
            dim = disc.dimension;
            B = cell(size(blocks,1),1);
            P = cell(size(blocks,1),1);
            for i = 1:size(blocks,1)
                M = cat(2,blocks{i,:});
                [B{i},P{i}] = resize(disc,M,dim-reduce(i),dim);
            end
        end
        
        function L = blockDiscretize(disc,block)
            if isa(block,'linBlock')
                L = block.stack( disc );
            elseif isa(block,'chebfun')
                L = disc.toValues(block);
                if ( block.isTransposed )
                    L = L.';
                end
            elseif isnumeric(block)
                L = block;
            else
                error('Unrecognized block type.')
            end
        end
       
    end
    
end