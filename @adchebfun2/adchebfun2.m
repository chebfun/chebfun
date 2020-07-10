classdef (InferiorClasses = {?chebfun2}) adchebfun2
    %ADCHEBFUN2   A class for supporting automatic differentiation in Chebfun2.
    %
    %   The ADCHEBFUN2 class allows chebop2 to compute variable coefficients of
    %   partial differential operators.
    %
    %   This class is not intended to be called directly by the end user.
    %
    %   V = ADCHEBFUN2(U), where U is a CHEBFUN2, returns the ADCHEBFUN2
    %   object V, which has a derivatives seeded as the identity
    %   operator on the domain of U.
    %
    %   V = ADCHEBFUN2(...), where the input is the same as would be
    %   passed to the CHEBFUN2 constructor, constructs a CHEBFUN2 from the
    %   input to the method. It then returns the ADCHEBFUN object V with
    %   the function part consisting of the CHEBFUN constructed, and the
    %   derivative part as the identity operator on the domain.
    %
    %   THIS CLASS HAS LIMITED FUNCTIONALITY AND IS DESIGNED SPECIFICALLY FOR
    %   CHEBOP2 REQUIREMENTS ONLY.
    %
    % See also ADCHEBFUN.
    
    % Copyright 2017 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    properties
        % FUNC: A CHEBFUN2 which corresponds to the function the
        % ADCHEBFUN2 represents.
        func
        
        % JACOBIAN: The Frechet derivative of the function the ADCHEBFUN2
        % represents, with respect to a selected basis variable. The basis
        % variable is selected at the start of computation with ADCHEBFUN2
        % objects, and has the identity operators as its Frechet derivative.
        jacobian
        
        % DOMAIN: Domain of the ADCHEBFUN2.
        domain
    end
    
    %% CLASS CONSTRUCTOR.
    methods
        
        % Main constructor. Convert a CHEBFUN2 to an ADCHEBFUN2.
        function V = adchebfun2(varargin)
            if ( nargin == 0 )
                % Return an empty ADCHEBFUN2 object.
                return
            elseif ( isa(varargin{1}, 'chebfun2') )
                V.func = varargin{:};  % Assign to the CHEBFUN2 field of V.
                V.jacobian = { 1 };
            else
                tempFunc = chebfun2(varargin{:});
                V = adchebfun2(tempFunc);
            end
        end
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        function fout = cos(fin)
            %COS   Cosine of an ADCHEBFUN2.
            
            fout = fin;
            fout.func = cos(fin.func);
            fout.jacobian = -sin(fin.func) * fin.jacobian;
            
        end
        
        function f = diff(f, varargin)
            %DIFF   DIFF for ADCHEBFUN2. See CHEBFUN2/DIFF for syntax.
            
            % Start by differentiating the CHEBFUN2 of the ADCHEBFUN2.
            f.func = diff(f.func, varargin{:});
            
            % Find out how many derivatives are required in each direction.
            if ( nargin == 1 ) % Default.
                nx = 0;
                ny = 1;
            elseif ( nargin == 2 )
                % Two arguments passed, so second argument must be the
                % order.
                order = varargin{1};
                if length(order) == 1 % DIFF in y is default.
                    nx = 0;
                    ny = order;
                elseif length(order) == 2
                    % Got a vector passed, first entry is number of derivatives in x
                    % direction, second entry is number of derivatives in y direction.
                    nx = order(1);
                    ny = order(2);
                else
                    error('CHEBFUN2:DIFF',...
                        'Undetermined direction of differentiation.');
                end
            else
                % Three arguments passed, second one is the order, third is
                % the dimension.
                order = varargin{1};
                dim = varargin{2};
                if ( dim == 1 )
                    nx = 0;
                    ny = order;
                else
                    nx = order;
                    ny = 0;
                end
            end
            
            % Update the derivative information.
            [mA, nA] = size(f.jacobian);
            
            % The output cell will have size (mA+yOrder)x(nA+xOrder)
            newCell = num2cell(zeros(mA + ny, nA + nx));
            
            % Differentiation in the x and y direction. This shifts derivative
            % information to the right and to the bottom, which amounts to
            % putting the input cell at the bottom right end side of the output
            % cell.
            newCell(end-mA+1:end,end-nA+1:end) = f.jacobian;
            
            % Assign newCell to the output
            f.jacobian = newCell;
            
        end
        
        function f = diffx(f, k)
            %DIFFX   DIFFX for ADCHEBFUN2.
            
            if ( nargin < 2 )
                k = 1;
            end
            f = diff(f, k, 2);
            
        end
        
        function f = diffy(f, k)
            %DIFFY   DIFFY for ADCHEBFUN2.
            
            if ( nargin < 2 )
                k = 1;
            end
            f = diff(f, k, 1);
            
        end
        
        function f = divergence( F )
            %DIVERGENCE   Divergence of an ADCHEBFUN2.
            
            f = diff(F(1), 1, 2) + diff(F(2), 1, 1);
            
        end
        
        function F = gradient(f)
            %GRADIENT   Gradient of an ADCHEBFUN2.
            
            fx = diff(f, 1, 2);   % DIFF in x-variable.
            fy = diff(f, 1, 1);   % DIFF in y-variable.
            F = [ fx, fy ];
            
        end
        
        function out = isempty(f)
            %ISEMPTY   True for an empty ADCHEBFUN2.
            
            out = isempty(f.func);
            
        end
        
        function L = lap(f) 
            % LAP   Shorthand for LAPLACIAN
            
            L = laplacian( f ); 
        end
        
        function L = laplacian(f)
            %LAPLACIAN   Laplacian of an ADCHEBFUN2.
            
            L = diff(f, 2, 2) + diff(f, 2, 1);
            
        end
        
        function f = minus(f, g)
            %-	  MINUS method for ADCHEBFUN2.
            
            % f-g =  f + (-g).
            f = plus(f, uminus(g));
            
        end
        
        function F = mrdivide(f, g)
            %/   Right scalar divide for ADCHEBFUN2.
            
            F = f;
            if ( isa(g, 'double') )
                F = F * ( 1 / g );
            else
                error('ADCHEBFUN2:mrdivide:ADchebfun2','Did you mean ./ ?');
            end
            
        end
        
        function h = mtimes(f, g)
            %*   ADCHEBFUN2 multiplication by scalar.
            
            if ( isempty(f) || isempty(g) )  
                % Return an empty ADCHEBFUN2.
                h = adchebfun2;
                return
            end
            
            if ( (isa(f, 'chebfun2') && isa(g, 'adchebfun2')) ...
                    || (isa(f, 'adchebfun2') && isa(g, 'chebfun2')) )
                % catch error: CHEBFUN2 * ADCHEBFUN2
                msg1 = 'Cannot use * between chebfun2 and adchebfun2.\n';
                msg2 = 'Did you mean .* ?';
                msg = strcat(msg1, msg2);
                error('ADCHEBFUN2:mtimes:InvalidOperation', msg);
            end
            
            if ( isa(f, 'double') ) % DOUBLE * ADCHEBFUN2.
                h = g;
                h.func = f * vertcat(g.func);          
                h.jacobian = cellfun(@mtimes, g.jacobian, ...
                repmat({f}, size(g.jacobian)), 'UniformOutput', false);
            elseif ( isa(g,'double') ) % ADCHEFUN2 * DOUBLE.
                h = mtimes(g, f);
            end
                
        end
        
        function varargout = plot(f)
            %PLOT   Plot an ADCHEBFUN2.
            
            if ( nargout > 0 )
                varargout = { plot(f.func) };
            else
                plot(f.func);
            end
            
        end
        
        function h = plus(f, g)
            %+   PLUS for ADCHEBFUN2.
            
            if ( isempty(f) || isempty(g) )
                % Return an empty ADCHEBFUN2.
                h = adchebfun2;
                return
            end
            
            if ( ~isa(f, 'adchebfun2') ) % First argument was not an ADCHEBFUN2.
                % Swap arguments.
                h = plus(g, f);
                return
            end
            
            if isa(g, 'adchebfun2') % ADCHEBFUN2 + ADCHEBFUN2.
                h = f;
                
                % Add the CHEBFUN2s.
                h.func = f.func + g.func;
                
                % Obtain the dimensions of the derivatives of the inputs.
                [mA, nA] = size(f.jacobian);
                [mB, nB] = size(g.jacobian);
                
                % The final derivative will have the dimensions corresponding to the
                % maximum dimensions of the input derivatives. We create a cell with
                % all zero entries of that dimensions, then add the input derivative
                % information to the top left of that matrix
                newA = num2cell(zeros(max(mA,mB), max(nA,nB)));
                newB = newA;
                
                % Replace entries of the cells of the correct size with information
                % from the inputs.
                newA(1:mA,1:nA) = f.jacobian;
                newB(1:mB,1:nB) = g.jacobian;
                
                % Add the two derivative cells together, and return.
                h = f;
                h.jacobian = cellfun(@plus, newA, newB, 'UniformOutput', false);
                
            elseif ( isa(g, 'chebfun2') || isa(g, 'double') ) 
                % ADCHEBFUN2 + CHEBFUN2 or DOUBLE.
                h = f;
                h.func = f.func + g;
            else
                error('ADCHEBFUN2:plus:type',...
                    'Cannot add these two objects together.');
            end
            
        end
        
        function fout = sin(fin)
            %SIN   Sine of an ADCHEBFUN2.
            
            fout = fin;
            fout.func = sin(fin.func);
            fout.jacobian = cos(fin.func) * fin.jacobian;
            
        end
        
        function f = times(f, g)
            %.*   ADCHEBFUN2 multiplication.
            
            if ( isempty(f) || isempty(g) )
                % Return an empty ADCHEBFUN2.
                f = chebfun2;   
                return
            end
            
            % ADCHEBFUN2 * DOUBLE or DOUBLE * ADCHEBFUN2: call MTIMES.
            if ( isa(f, 'double') || isa(g, 'double') )
                f = mtimes(f, g);
                
            % ADCHEBFUN2 * CHEBFUN2 or CHEBFUN2 * ADCHEBFUN2 does not need complicated
            % derivative computation, as a CHEBFUN2 does not have any derivative
            % information, and can just be treated like a constant.
            elseif ( isa(f, 'adchebfun2') && isa(g, 'chebfun2') )
                f = times(g, f);   
            elseif ( isa(f, 'chebfun2') && isa(g, 'adchebfun2') )
                g.func = f .* g.func;
                g.jacobian = cellfun(@mtimes, g.jacobian, ...
                    repmat({f}, size(g.jacobian)), 'UniformOutput', false);
                f = g;   
            else
                error('ADCHEBFUN2:times',...
                    'Can only multiply an adchebfun2 by a scalar or a chebfun2.');
            end
            
        end
        
        function f = uplus(f)
            %+   Unary plus for ADCHEBFUN2.
            
        end
        
        function f = uminus(f)
            %-   Unary minus for ADCHEBFUN2.
            
            f.func = uminus(f.func);
            f.jacobian = cellfun(@uminus, f.jacobian, 'UniformOutput', false);
            
        end
        
    end
    
end
