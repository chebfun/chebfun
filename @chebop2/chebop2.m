classdef chebop2
%CHEBOP2   CHEBOP2 class for representing linear partial differential operators.
%
% Class used to solve linear PDEs defined on rectangular domains that have unique and
% globally smooth solutions.
%
% N = CHEBOP2(@(u) op(u)) constructs an operator N representing the operator
% given in @(u)op(u) acting on functions of two variables on [-1,1] by [-1,1].
%
% N = CHEBOP2(@(u) op(u), [a b c d]) constructs an operator N acting on
% functions of two variables defined on [a,b] by [c,d].
%
% N = CHEBOP2(@(x,y,u) op(x,y,u),...) constructs a variable coefficient PDE
% operator. If a partial differential operator is constant coefficient, we
% recommend the @(u) notation rather than @(x,y,u) as it is more efficient. 
%
% Boundary conditions are imposed via the syntax N.lbc, N.rbc, N.ubc, and N.dbc.
%
% Example 1: (Poisson with Dirichlet conditions):
%    N = chebop2(@(u) diff(u,2,1) + diff(u,2,2));
%    N.bc = 0;
%    u = N \ 1;
% 
% Example 2: (Helmholtz equation with gravity)
%    N = chebop2(@(x,y,u) laplacian(u) - 10*y.^2.*u, [-1 1 -3 0]); 
%    N.bc = 1; 
%    u = N \ 0; 
%
% Example 3: (Klein-Gordon equation) 
%    N = chebop2(@(u) diff(u,2,1) - diff(u,2,2) + 5*u,[-1 1 0 3]); 
%    N.lbc = 0; N.rbc = 0; 
%    N.dbc = @(x,u) [u - exp(-30*x.^2) ; diff(u)];
%    u = N \ 0; 
% 
% For further details about the PDE solver, see:
%
% A. Townsend and S. Olver, The automatic solution of partial differential
% equations using a global spectral method, in preparation, 2014.
%
% Warning: This PDE solver is an experimental new feature. It has not been
% publicly advertised.  Chebop2 cannot do nonlinear problems as more 
% algorithmic advances are needed. 

% Copyright 2016 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %% PROPERTIES.
    properties ( GetAccess = 'public', SetAccess = 'public' )
        
        domain = [];  % Domain of the operator.
        op = [];      % The operator.
        ubc = [];     % Up boundary condition(s).
        lbc = [];     % Left boundary condition(s).
        rbc = [];     % Right boundary condition(s).
        dbc = [];     % Down boundary condition(s).
        dim = [];     % Size of the system (number of eqns).
        scale = [];   % Relative solution scale.
        coeffs = [];  % Matrix storing constant coefficients.
        rhs = [];     % Righthand side, if given by user. 
        xorder = 0;   % Diff order in the x-variable.
        yorder = 0;   % Diff order in the y-variable.
        U             %
        S             % Low rank form of the partial differential operator.
        V             %
        
    end
    
    %% CONSTRUCTOR.
    methods
        
        function N = chebop2(varargin)
            % CHEBOP2 CONSTRUCTOR.
            
            % Get CHEBFUN2 preferences.
            pref = chebfunpref();
            tol = pref.cheb2Prefs.chebfun2eps;
            
            % If empty input arguments then return an empty CHEBOP2 object.
            if ( isempty(varargin) )
                return
            end
            
            % What domain is the operator defined on?
            if ( numel(varargin) > 1 )
                ends = varargin{2}; % Second argument should be a domain.
                if ( length(ends) == 4 )
                    % Valid domain?
                    if ( diff(ends(1:2)) > 0 && diff(ends(3:4)) > 0 )
                        dom = ends;
                    else
                        error('CHEBFUN:CHEBOP2:chebop2:emptyDomain', ...
                            'Empty domain.');
                    end
                else
                    error('CHEBFUN:CHEBOP2:chebop2:badDomain',...
                        'Argument should be a domain given by four doubles.')
                end
            else
                if ( isa(varargin{1}, 'function_handle') )
                    % Pick the default domain.
                    rect1 = [-1, 1];
                    rect2 = [-1, 1];
                    dom = [rect1, rect2];
                elseif ( isa( varargin{1}, 'double') )
                    % Set up identity operator on the domain.
                    N = chebop2(@(u) u, varargin{1});
                    return
                else
                    error('CHEBFUN:CHEBOP2:chebop2:badArg',...
                        'First argument is not an operator or domain.')
                end
            end
            
            % First argument in the constructor is the operator. If the
            % operator is univariate then it's a constant coefficient PDE,
            % otherwise assume it is a variable coefficient.
            if ( isa(varargin{1}, 'function_handle') )
                fh = varargin{1};
                
                if ( nargin(fh) == 1 )  % The PDE has constant coefficients.
                    % Trust that the user has formed the CHEBFUN2 objects
                    % outside of CHEBOP2.
                    
                    % Extract out rhs: 
                    x = chebfun2(@(x,y) x, dom); 
                    RHS = fh(0*x);
                    % In this case, the RHS must be a constant chebfun2.
                    if ( norm( RHS ) > 0 ) 
                        fh = @(u) fh(u) - RHS; 
                        N.rhs = -RHS;             % store for later.
                    else
                        N.rhs = 0; 
                    end

                    u = adchebfun2(chebfun2(@(x,y) x.*y, dom));
                    v = fh(u);
                    % If the PDO has constant coefficients then convert to
                    % double:
                    try
                        A = cell2mat(v.jacobian).';
                    catch
                        % PDO has variable coefficients, keep them in a
                        % cell array:
                        A = v.jacobian;
                    end
                    
                elseif ( nargin(fh) == 2 )
                    error('CHEBFUN:CHEBOP2:chebop2:badOp1', ...
                        'Did you intend to have @(x,y,u)?')
                elseif ( nargin(fh) == 3 )
                    % The coefficients of the PDE are now variable
                    % coefficient.
                    
                    % Setup a chebfun2 on the right domain
                    u = adchebfun2(chebfun2(@(x,y) x.*y, dom));
                    x = chebfun2(@(x,y) x, dom);
                    y = chebfun2(@(x,y) y, dom);
                    
                    % Extract out rhs: 
                    RHS = fh(x, y, 0*x); 
                    % The RHS is a chebfun2. Check if it is the zero. Do
                    % not change function handle unless we have to because 
                    % we want the user to be able to see the PDO 
                    % untouched if possible. 
                    if ( norm( RHS ) > 0 ) 
                        fh = @(x, y, u) fh(x, y, u) - RHS; 
                        N.rhs = -RHS;             % store for later.
                    else
                        N.rhs = 0; 
                    end
                    
                    % Apply it to the operator.
                    v = fh(x, y, u);
                    A = v.jacobian;  % Cell array of variable coefficients.
                    
                    % If we have a variable coefficient PDO, then compute the 
                    % separable representation immediately. We need it now.
                    [cellU, matS, cellV] = chebop2.separableFormat( A,...
                                                    size(A,2), size(A,1), dom );
                    N.U = cellU;
                    N.S = matS;
                    N.V = cellV;
                else
                    error('CHEBFUN:CHEBOP2:chebop2:badOp2',...
                        'Operator should be @(u) or @(x,y,u).')
                end
                
            else
                error('CHEBFUN:CHEBOP2:chebop2:badOp3',...
                    'First argument should be an operator')
            end
            
            % Often the coefficients are obtained with small rounding errors
            % and it is important to remove the very small non-zero ones to
            % have rank(A) correct.
            if ( iscell(A) )
                for jj = size(A, 1)
                    for kk = size(A, 2)
                        if ( isa(A{jj,kk}, 'double') && abs(A{jj,kk}) < 10*tol )
                            A{jj,kk} = 0;
                        end
                    end
                end
            else
                A(abs(A) < 10*tol) = 0;
            end
            
            % Construct CHEBOP2 object. The boundary conditions will be
            % given later.
            N.domain = dom;
            N.op = fh;
            N.coeffs = A;
            
            % Calculate xorder and yorder of PDE.
            % Find the differential order of the PDE operator.
            if ( iscell(A) )
                xdifforder = size(A, 2) - 1;
                ydifforder = size(A, 1) - 1;
            elseif ( min(size(A)) > 1 )
                xdifforder = find(sum(abs(A), 2) > 100*tol, 1, 'last') - 1;
                ydifforder = find(sum(abs(A)) > 100*tol, 1, 'last' ) - 1;
            else
                if ( size(A, 1) == 1 )
                    ydifforder = length(A) - 1;
                    xdifforder = 0;
                else
                    xdifforder = length(A) - 1;
                    ydifforder = 0;
                end
            end
            N.xorder = xdifforder;
            N.yorder = ydifforder;
            
            % Issue a warning to the user for the first CHEBOP2:
            warning('CHEBFUN:CHEBOP2:chebop2:experimental',...
                ['CHEBOP2 is a new experimental feature.\n'...
                'It has not been tested to the same extent as other'...
                'parts of the software.']);
            % Turn it off:
            warning('off', 'CHEBFUN:CHEBOP2:chebop2:experimental');
            
        end
        
    end
    
    %% STATIC HIDDEN METHODS.
    methods ( Static = true, Hidden = true )
        
        % Matrix equation solver: AXB^T + CXD^T = E. xsplit, ysplit = 1 if
        % the even and odd modes (coefficients) decouple.
        X = bartelsStewart(A, B, C, D, E, xsplit, ysplit);
        
        % Use automatic differentiation to pull out the coeffs of the
        % operator:
        deriv = chebfun2deriv(op);
        
        % This is used to discretize the linear constrains:
        [bcrow, bcvalue] = constructBC(bcArg, bcpos,...
            een, bcn, dom, scl, order);
        
        % This is used to discretize the PDE:
        [CC, rhs, bb, gg, Px, Py, xsplit, ysplit] =...
            discretize(N, f, m, n, flag);
        
        % Convert all the different user inputs for bc into uniform format:
        bc = createBC(bcArg, ends);
        
        % Method for deciding how to solve the matrix equation:
        X = denseSolve(N, f, m, n);
        
        % Compute the separable representation of a PDO: 
        [cellU, S, cellV] = separableFormat(A, xorder, yorder, dom);
        
        % Remove trailing coefficients.
        a = truncate(a, tol);
        
    end
    
end
