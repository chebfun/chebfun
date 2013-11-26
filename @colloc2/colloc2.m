classdef colloc2 < chebDiscretization
     
    methods
        function disc = colloc2(source,dimension,domain)
            if isempty(source)
                return
            end
            disc.source = source; 
            disc.domain = source.domain;

            if nargin > 1
                disc.dimension = dimension;
                if nargin > 2
                    disc.domain = domain;
                end
            end
        end
         
    end
    
    methods
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
                                       
        function f = toFunction(disc,values)
            if ( disc.numIntervals > 1 )
                values = mat2cell(values,disc.dimension);
            end
            f = chebfun(values, disc.domain);
        end
        
    end
    
    methods
        
        function varargout = matrix(dsc,dimension,domain)
            if nargin > 1
                dsc.dimension = dimension;
                if nargin > 2
                    dsc.domain = domain;
                end
            end
            L = dsc.source;
            if isa(L,'chebmatrix')
                A = cellfun(@(x) blockMatrix(dsc,x),L.blocks,'uniform',false);
                out{1} = A;
                if isa(L,'linop')
                    [out{1:3}] = toLinop(dsc,A);
                end
                m = max(1,nargout);
                varargout(1:m) = out(1:m);
            else
                [varargout{1:nargout}] = blockMatrix(dsc);
            end
        end
            
        function [A,P,B] = toLinop(disc,blocks)
            L = disc.source;
            
            [rows,P] = disc.reproject(blocks);
            P = blkdiag(P{:});

            B = [];
            if ~isempty(L.constraint)
                disc.source = L.constraint.operator;
                constr = matrix(disc);
                B = [ cell2mat(constr); B ];
            end
            if ~isempty(L.continuity)
                disc.source = L.continuity.operator;
                constr = matrix(disc);
                B = [ cell2mat(constr); B ];
            end
            A = cell2mat(rows);
            A = [ B; A ];

        end
        
        function A = blockMatrix(disc,item)
            if ( nargin < 2 )
                item = disc.source;
            end
            if isa(item,'linBlock')
                disc.source = item;
                A = disc.source.stack( disc );
            elseif isa(item,'chebfun')
                A = disc.toValues(item);
                if ( item.isTransposed )
                    A = A.';
                end
            elseif isnumeric(item)
                A = item;
            else
                error('Unrecognized block type.')
            end
        end
        
        
        function [v,disc] = mldivide(disc,A,b)
            s = 1./ max(1, max(abs(A),[],2) );
            A = bsxfun(@times,s,A);
            v = A \ (s.*b);
        end
        
        function b = rhs(disc,f)
            if isempty(disc.dimension)
                error('Discretization dimension not given.')
            end
            % NONONO
            if isa(f,'chebfun'), f = chebmatrix({f}); end
            row = cellfun(@(x) blockMatrix(disc,x),f.blocks,'uniform',false);
%             if ( ~iscell(row) )
%                 row = {row};
%             end
            row = disc.reproject(row);
            
            b = cell2mat(row);
            L = disc.source;
            if ~isempty(L.constraint)
                b = [ L.constraint.values; b ];
            end
            if ~isempty(L.continuity)
                b = [ L.continuity.values; b ];
            end
        end
        
    end
    
    methods ( Access = private )
        function [B,P] = reproject(disc,blocks)
            reduce = disc.source.sizeReduction;
            dim = disc.dimension;
            B = cell(size(blocks,1),1);
            P = cell(size(blocks,1),1);
            for i = 1:size(blocks,1)
                M = cat(2,blocks{i,:});
                [B{i},P{i}] = resize(disc,M,dim-reduce(i),dim);
            end
        end
        
%         function L = blockDiscretize(disc,block)
%             if isa(block,'linBlock')
%                 L = block.stack( disc );
%             elseif isa(block,'chebfun')
%                 L = disc.toValues(block);
%                 if ( block.isTransposed )
%                     L = L.';
%                 end
%             elseif isnumeric(block)
%                 L = block;
%             else
%                 error('Unrecognized block type.')
%             end
%         end
       
    end
    
end