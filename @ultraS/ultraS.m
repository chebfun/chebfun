classdef ultraS < chebDiscretization
    properties
        coeffs
        outputSpace
    end
    
    methods
        function disc = ultraS(source,dimension,domain)
            %ULTRAS constructor.
            
            if ( nargin == 0 || isempty(source) )
                % Construct an empty ULTRAS.
                return
            end
            
            if ( nargin > 1 )
                disc.dimension = dimension;
                if ( nargin > 2 )
                    disc.domain = domain;
                end
            end
            
%             % Decide which kind of object we're discretizing:
%             if ( isa(source, 'linop') )
%                 disc.linop = source;
%             else
                disc.source = source;
%            end
            
            % Obtain the coeffs and output psace required for this source:
            disc.coeffs = disc.getCoeffs(source);
            disc.outputSpace = disc.getOutputSpace(source);
            
        end
        
        function L = blockDiscretize(disc, block)
            if (isa(block, 'operatorBlock') )
                if ( ~isempty(disc.coeffs) )
                    L = ultraS.quasi2USdiffmat(disc.coeffs, disc.domain, ...
                        disc.dimension, disc.outputSpace);
                else
                    error
                end
            elseif ( isa(block, 'functionalBlock') )
                dim = disc.dimension;
                dom = disc.domain;
                collocDisc = colloc2(block,dim,dom);
                L = matrix(collocDisc);
                cumsumDim = [0, cumsum(dim)];
                numInts = numel(dom) - 1;
                tmp = cell(1, numInts);
                for k = 1:numInts
                    Lk = L(cumsumDim(k) + (1:dim(k)));
                    tmp{k} = flipud(chebtech2.coeffs2vals(Lk.')).';
                end
                L = cell2mat(tmp);
            elseif ( isa(block, 'chebfun') )
                L = toValues(disc, block);
                if ( block.isTransposed )
                    error % TODO: ?
                    L = L.';
                end
            elseif ( isnumeric(block) )
                L = block;
            else
                error('Unrecognized block type.')
            end
        end
        
        function S = convert( A, K1, K2 )
            %CONVERT(A, K1, K2), convert C^(K1) to C^(K2)
            d = A.domain;
            n = A.dimension;
            numIntervals = length(d) - 1;
            % Find the diagonal blocks.
            blocks = cell(numIntervals);
            for k = 1:numIntervals
                blocks{k} = A.convertmat(n(k), K1, K2);
            end
            % Assemble.
            S = blkdiag(blocks{:});
        end
        
        function D = diff(A, m)
            d = A.domain;
            n = A.dimension;
            if ( m == 0 )
                D = speye(sum(n));
            else
                numIntervals = length(d) - 1;
                % Find the diagonal blocks.
                blocks = cell(numIntervals);
                for k = 1:numIntervals
                    len = d(k+1) - d(k);
                    blocks{k} = A.diffmat(n(k), m) * (2/len)^m;
                end
                % Assemble.
                D = blkdiag(blocks{:});
            end
        end
        
        function L = matrix(disc,dimension,domain)
            % TODO: error checking on inputs
            if ( nargin > 1 )
                disc.dimension = dimension;
                if ( nargin > 2 )
                    disc.domain = domain;
                end
            end
            
            A = disc.source;
%            validate(disc);
            if ( isa(A, 'chebmatrix') )
                %diffOrder = getDiffOrder(A);  % TODO: not used?
                c = disc.coeffs;
                outputSpaces = disc.outputSpace;
                L = cell(size(A));
                for j = 1:size(A, 1)
                    disc.outputSpace = outputSpaces(j);
                    for k = 1:size(A, 2)
                        disc.coeffs = c{j,k};
                        L{j,k} = blockDiscretize(disc, A.blocks{j,k});
                    end
                end
                if ( isa(A,'linop') )
                    L = toLinop(disc,L);
                end
            else
                disc.coeffs = disc.coeffs{1};
                L = blockDiscretize(disc, A);
            end
        end
        
        function A = toLinop(disc,blocks)
            L = disc.source;
            rows = disc.reproject(blocks);
            A = cell2mat(rows);
            dim = disc.dimension;
            dom = disc.domain;
            if ~isempty(L.constraint)
                disc2 = ultraS(L.constraint.operator,dim,dom);
                constr = matrix(disc2);
                A = [ cell2mat(constr); A ];
            end
            if ~isempty(L.continuity)
                disc2 = ultraS(L.continuity.operator,dim,dom);
                constr = matrix(disc2);
                A = [ cell2mat(constr); A ];
            end
        end
        
        function B = reproject(disc, blocks)
            reduce = disc.source.sizeReduction;
            dim = disc.dimension;
            B = cell(size(blocks,1),1);
            for i = 1:size(blocks,1)
                M = cat(2, blocks{i,:});
                B{i} = resize(disc, M, dim - reduce(i));
            end
        end
        
        function [v,disc] = mldivide(disc,A,b)
            v = A\b;
        end
 
        function A = resize(disc, A, m)
            dom = disc.domain;
            n = disc.dimension;
            % chop off some rows and columns
            v = [];
            nn = cumsum([0 n]);
            for k = 1:numel(dom)-1
                v = [v m(k) + nn(k) + (1:(n(k)-m(k)))];
            end
            A(v.', :) = [];
        end
        
        function b = rhs(disc, f)
            if isempty(disc.dimension)
                error('Discretization dimension not given.')
            end
            fdisc = ultraS(f,disc.dimension,disc.domain);
            row = matrix(fdisc);
            if ~iscell(row)
                row = {row};
            end
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
        
        function f = toFunction(disc, values)
            dom = disc.domain;
            v = mat2cell(full(values), disc.dimension, 1);
            funs = cell(numel(v),1);
            for k = 1:numel(v)
                ct = chebtech2({[], flipud(v{k})});
                funs{k} = bndfun(ct, dom(k:k+1));
            end
            f = chebfun(funs);
        end
        
        function fx = toValues(disc, f)
            dom = disc.domain;
            numInts = numel(dom) - 1;
            dim = disc.dimension;
            f = restrict(f, dom);
            c = cell(numInts, 1);
            for k = 1:numInts
                c{k} = flipud(chebpoly(f, k, dim(k)).');
            end
            c = cell2mat(c);
            S = convert(disc, 0, disc.outputSpace);
            fx = S*c;
            
        end
        
    end
    
    methods (Static)
        
        function S = convertmat( n, K1, K2 )
            %CONVERTMAT(A, K1, K2), convert C^(K1) to C^(K2)
            S = speye(n);
            for s = K1:K2
                S = spconvert(n, s) * S;
            end
        end
        
        function D = diffmat( n, m )
            %DIFFMAT(N, K, N), computes the kth order US derivative matrix
            if ( m > 0 )
                D = spdiags((0:n-1)', 1, n, n);
                for s = 1:m-1
                    D = spdiags(2*s*ones(n, 1), 1, n, n)*D;
                end
            else
                D = speye(n);
            end
        end
        
        function f_coeffs = discretizeFunction(f, dim, dom)
            if ( nargin < 3 )
                dom = f.domain;
            end
            f_coeffs = [];
            f = restrict(f, dom);
            for k = 1:numel(dom)-1
                dimk = dim(k);
                tmp = flipud(get(f.funs{k}, 'coeffs'));
                n = length(tmp);
                % prolong/truncate.
                if ( n > dimk )
                    tmp = tmp(1:dimk);
                else
                    tmp = [tmp ; zeros(dimk - n, 1)];
                end
                f_coeffs = [f_coeffs ; tmp];
            end
        end
        
        
        function c = getCoeffs(source)
            if ( isa(source, 'chebmatrix') )
                c = cell(size(source));
                for k = 1:numel(c)
                    try
                        c{k} = blockCoeff(source.blocks{k}).coeffs;
                    catch
                        c{k} = [];
                    end
                end
            else
                try
                    c = {blockCoeff(source).coeffs};
                catch
                    c = {};
                end
            end
        end
        
        function outputSpace = getOutputSpace(source)
            if ( isa(source, 'chebmatrix') )
                diffOrders = getDiffOrder(source);
            else
                diffOrders = source.diffOrder;
            end
            outputSpace = max(max(diffOrders, [], 2) - 1, -1);
        end
        
        function f = makeChebfun(u, dom)
            funs = cell(numel(u), 1);
            for k = 1:numel(u)
                ct = chebtech2({[], flipud(u{k})});
                funs{k} = bndfun(ct, dom(k:k+1));
            end
            f = chebfun(funs);
            u = chebmatrix({f});
        end
        
        L = quasi2USdiffmat(L, dom, dim, outputSpace)
        
        function [isDone, epsLevel] = testConvergence(v)
            % Test convergence on a single interval:
            v = full(v);
            f = chebtech2({[], flipud(v)});
            [isDone, epsLevel] = strictCheck(f);
        end
    end
end
