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
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    properties
        chebfun2   % An adchebfun2 has a chebfun2
        der        % Derivative information stored here
    end
    
    %% CLASS CONSTRUCTOR:
    methods
        % Main constructor. Convert a chebfun2 to ADchebfun2
        function g = adchebfun2 ( varargin )
            if( nargin == 0 )
                % return an empty chebfun2 object.
            elseif isa(varargin{1},'chebfun2')
                g.chebfun2 = varargin{:};  % Assign to the chebfun2 field of g.
                g.der = chebfun2der(1,g.chebfun2.domain);
            else
                cheb2temp = chebfun2(varargin{:});
                g = adchebfun2(cheb2temp);
            end
        end
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        function fout = cos(fin)
            %COS Cosine of ADCHEBFUN2.
            fout = fin;
            fout.chebfun2 = cos(fin.chebfun2);
            fout.der = (-sin(fin.chebfun2))*fin.der;
        end
        
        function f = diff(f, varargin)
            % DIFF for ADCHEBFUN2. See chebfun2/diff for syntax.
            
            % Start by differentiating the chebfun2 of the ADchebfun
            f.chebfun2 = diff(f.chebfun2, varargin{:});
            
            % Find out how many derivatives are required in each direction
            if ( nargin == 1 ) % defaults.
                nx = 0;
                ny = 1;
            elseif ( nargin == 2 )
                % Two arguments passed, so second argument must be the order
                order = varargin{1};
                if length(order) == 1 % diff in y is default.
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
                % Three arguments passed, second one is the order, third is the dimension
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
            
            % Update the derivative information
            f.der = diff(f.der, nx, ny);
        end
        
        function f = diffx(f, k)
            % DIFFX for ADCHEBFUN2
            if ( nargin < 2 )
                k = 1;
            end
            f = diff( f, k, 2 );
        end
        
        function f = diffy(f, k)
            % DIFFY for ADCHEBFUN2
            if ( nargin < 2 )
                k = 1;
            end
            f = diff(f, k, 1);
        end
        
        function f = divergence( F )
            %DIVERGENCE   Divergence of ADCHEBFUN2.
            f = diff( F(1), 1, 2 ) + diff( F(2), 1, 1 );
        end
        
        function F = gradient( f )
            %GRADIENT  Numerical gradient of ADCHEBFUN2.
            fx = diff(f, 1, 2);   % diff in x-variable
            fy = diff(f, 1, 1);   % diff in y-variable
            F = [ fx fy ];
        end
        
        function L = laplacian( f )
            %LAPLACIAN   Laplacian of ADCHEBFUN2.
            
            % laplacian(f) = f_xx + f_yy:
            L = diff(f, 2, 2) + diff(f, 2, 1);
        end
        
        function f = minus(f, g)
            % -	  Minus method for ADCHEBFUN2
            
            % f-g =  f + (-g)
            f = plus(f,uminus(g));
        end
        
        function F = mrdivide(f, g)
            % /   Right scalar divide for ADCHEBFUN2.
            
            F = f;
            if (isa(g,'double'))
                F = F * ( 1 / g );
            else
                error('ADCHEBFUN2:mrdivide:ADchebfun2','Did you mean ./ ?');
            end
        end
        
        function h = mtimes(f,g)
            %*	ADchebfun2 multiplication.
            
            if ( isempty(f) || isempty(g) )  % just return an empty chebfun2
                h = chebfun2;
                return
            end
            
            if ( isa(f,'double') )
                
                % double * chebfun2
                h = g;
                
                % Multiply the chebfun2
                h.chebfun2 = f * ( vertcat( g.chebfun2  ));
                
                % Multiply the derivative
                h.der = ( g.der ) * f;
                
            elseif ( isa(g,'double') )
                % chebfun2 * double
                h = f;
                
                % Multiply the chebfun2
                h.chebfun2 = g*(f.chebfun2);
                
                % Multiply the derivative
                h.der = (f.der)*g;
            end
        end
        
        function varargout = plot(f)
            % Plot an adchebfun2
            if ( nargout )
                varargout = plot(f.chebfun2);
            else
                plot(f.chebfun2);
            end
        end
        
        function h = plus ( f, g )
            %+	  Plus for ADCHEBFUN2
            
            if ( isempty(f) || isempty(g) )  % check for empty ADchebfun2
                h = ADchebfun2;  % just return an empty ADchebfun2.
                return
            end
            
            if ( ~isa(f, 'adchebfun2') ) % First argument was not an ADchebfun2
                % Swap arguments.
                h = plus(g, f);
                return
            end
            
            if isa(g, 'adchebfun2')     % ADCHEBFUN2 + ADCHEBFUN2
                h = f;
                
                % Add the chebfun2s
                h.chebfun2 = f.chebfun2 + g.chebfun2;
                
                % Add the der fields
                h.der = f.der + g.der;
                
            elseif ( isa(g, 'chebfun2') || isa(g, 'double') ) % ADCHEBFUN2 + DOUBLE/SCALAR
                f.chebfun2 = f.chebfun2 + g;
            else
                error('ADCHEBFUN2:plus:type',...
                'Cannot add these two objects together');
            end
            
        end
        
        function f = sin( f )
            % SIN   Sine of an ADCHEBFUN2.
            fout = fin;
            fout.chebfun2 = sin(fin.chebfun2);
            fout.der = (cos(fin.chebfun2))*fin.der;
        end
        
        function f = times(f, g)
            % .*   ADChebfun2 multiplication.
            
            if ( isempty( f ) || isempty( g ) )
                f = chebfun2;   % just return an empty chebfun2.
                return
            end
            
            % ADCHEBFUN2*DOUBLE or DOUBLE*ADCHEBFUN2. Just call mtimes.
            if ( isa(f, 'double') || isa(g, 'double') )
                f = mtimes(f,g);
                
                % ADCHEBFUN2*CHEBFUN2 or CHEBFUN2*ADCHEBFUN2 does not need complicated
                % derivative computation, as a CHEBFUN2 does not have any derivative
                % information, and can just be treated like a constant
            elseif ( isa(f, 'adchebfun2') && isa(g, 'chebfun2') )
                f.chebfun2 = f.chebfun2.*g;
                f.der = f.der*g;
            elseif ( isa(f, 'chebfun2') && isa(g, 'adchebfun2') )
                g.chebfun2 = f .* g.chebfun2;
                g.der = ( g.der ) * f;
                f = g;
            elseif ( isa(f,'chebfun2') && isa(g,'chebfun2v') )
                %% chebfun2 * chebfun2v
                
                % This functionality may be taken out of a release.
                g.xcheb = f.*g.xcheb;
                g.ycheb = f.*g.ycheb;
                if ~isempty(g.zcheb)
                    g.zcheb = f.*g.zcheb;
                end
                f = g;
            else
                % We had a chebfun2.*unknown, so complain.
                error('Chebfun2:times',...
                    'Can only do chebfun2 times scalar or chebfun2.');
            end
        end
        
        function f = uplus( f ) 
            % +  Unary plus for ADCHEBFUN2.
        end
        
        function f = uminus( f )
            % -	  Unary minus for ADCHEBFUN2.
            f.chebfun2 = uminus( f.chebfun2 );
            f.der = -f.der;
        end
    end
end
