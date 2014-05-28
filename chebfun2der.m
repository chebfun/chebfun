classdef  (InferiorClasses = {?chebfun2}) chebfun2der
    %CHEBFUN2DER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        derCell;
        domain;
    end
    
    methods
        % Main constructor. Convert scalar input to a cell
        function der = chebfun2der (value,dom)
            if( nargin == 0 )
                % return an empty chebfun2der.
            else
                der.derCell = {value};  % Assign to the derCell field of der.
                der.domain = dom;
            end
        end
        
        function C = plus(A, B)
            
            % Obtain the dimensions of the derivatives of the inputs.
            [mA, nA] = size(A.derCell);
            [mB, nB] = size(B.derCell);
            
            % The final derivative will have the dimensions corresponding to the
            % maximum dimensions of the input derivatives. We create a cell with
            % all zero entries of that dimensions, then add the input derivative
            % information to the top left of that matrix
            newA = num2cell(zeros(max(mA,mB),max(nA,nB)));
            newB = newA;
            
            % Replace entries of the cells of the correct size with information
            % from the inputs.
            newA(1:mA,1:nA) = A.derCell;
            newB(1:mB,1:nB) = B.derCell;
            
            % Add the two derivative cells together, and return
            C = A;
            C.derCell = cellfun(@plus, newA,newB,'UniformOutput',false);
        end
        
        function A = uminus(A)
            % uminus via cell-fun
            A.derCell = cellfun(@uminus,A.derCell,'UniformOutput',false);
        end
        
        function A = diffX(A , xOrder)
            A = diff(A, xOrder, 0);
        end
        
        function A = diffY(A ,yOrder)
            A = diff(A, 0, yOrder);
        end
        
        function A = diff(A, xOrder, yOrder)
            % Obtain the size of the input derivative
            [mA, nA] = size(A.derCell);
            
            % The output cell will have size (mA+yOrder)x(nA+xOrder)
            newCell = num2cell(zeros(mA + yOrder, nA + xOrder));
            
            % Differentiation in the x and y direction. This shifts derivative
            % information to the right and to the bottom, which amounts to
            % putting the input cell at the bottom right end side of the output
            % cell.
            newCell(end-nA+1,end-mA+1:end) = A.derCell;
            
            % Assign newCell to the output
            A.derCell = newCell;
        end
        
        function A = mtimes(A,s)
            % Can assume A is a chebfun2der and s is a scalar, since that's the
            % only way those arguments are going to be passed to this method
            
            % Swap arguments if first argument is not a chebfun2der double
            if ~isa(A,'chebfun2der')
                A = mtimes(s,A);
                return
            end
            
            % Can assume A is now a chebfun2der, s is either a scalar, or a
            % chebfun2
            
            if isa(s,'double') || isa(s,'chebfun2')
                %                 A.derCell = cellfun(@mtimes,A.derCell, ...
                %                     num2cell(s*ones(size(A.derCell))), 'UniformOutput', false);
                %             else isa(s,'chebfun2')
                A.derCell = cellfun(@mtimes,A.derCell, ...
                    repmat({s},size(A.derCell)), 'UniformOutput', false);
            end
        end
        
        function h = plot(A)
            Acell = A.derCell;
            Adom = A.domain;
            [m,n]=size(Acell);
            kk = 1;
            for j = 1:m
                for k = 1:n
                    f = Acell{j,k};
                    subplot(m,n,kk)
                    if isa(f,'double')
                        xFill = [Adom(1) Adom(1) Adom(2) Adom(2)];
                        yFill = [Adom(3) Adom(4) Adom(4) Adom(3)];
                        if f == 0
                           fill(xFill,yFill,[255 105 105]./256); 
                        else
                            fill(xFill,yFill,[128 208 242]./256);
                        end
                    else
                        contour(f)
                    end
                    strTitle = 'u';
                    if k > 1
                        strTitle = [strTitle, '_{' repmat('x',[1 k-1]) '}'];
                    end
                    if j > 1
                        strTitle = [strTitle, '_{' repmat('y',[1 j-1]) '}'];
                    end
                    title(strTitle,'FontSize',14)
                    axis tight
                    kk = kk + 1;
                end
            end
        end
        
    end
    
end