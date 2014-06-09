classdef chebop2
%CHEBOP2  CHEBOP2 class for representing partial differential equations
%
% Class used to solve PDEs defined on rectangular domains that have
% unique and globally smooth solutions.
%
% N = CHEBOP2(@(u) op(u)) constructs an operator N representing the
% operator given in @(u)op(u) acting on functions of two variables on
% [-1,1] by [-1,1].
%
% N = CHEBOP2(@(u) op(u), [a b c d]) constructs an operator N acting on
% functions of two variables defined on [a,b] by [c,d].
%
% N = CHEBOP2(@(x,y,u) op(x,y,u),...) constructs a variable coefficient PDE
% operator.
%
% Boundary conditions are imposed via the syntax N.lbc, N.rbc, N.ubc, and
% N.dbc. For example to solve Poisson with Dirichlet conditions try:
%
% Example:
%    N = chebop2(@(u) diff(u,2,1) + diff(u,2,2));
%    N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0;
%    u = N \ 1;
%
% For further details about the PDE solver, see: 
% 
% A. Townsend and S. Olver, The automatic solution of partial differential
% equations using a global spectral method, in preparation, 2014.
% 
% Warning: This PDE solver is an experimental new feature. It has not been
% publicly advertised.  
        
% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        opshow = [];
        ubc = [];       % Up boundary condition(s)
        ubcshow=[];
        lbc = [];       % Left boundary condition(s)
        lbcshow=[];
        rbc = [];       % Right boundary condition(s)
        rbcshow=[];
        dbc = [];       % Down boundary condition(s)
        dbcshow=[];
        dim = [];       % Size of the system (number of eqns)
        scale = [];     % Relative solution scale
        coeffs=[];      % matrix storing constant coefficients
        xorder = 0;
        yorder = 0;
        U
        S               %explict low rank form.
        V
    end
    
    methods
        
        function N = chebop2( varargin )
            N = constructor( N, varargin{:} );
        end
        
    end
    
    methods ( Static = true, Hidden = true )
        % Matrix equation solver: AXB^T + CXD^T = E. xsplit, ysplit = 1 if
        % the even and odd modes (coefficients) decouple.
        X = BartelsStewart(A, B, C, D, E, xsplit, ysplit);
        
        % Use automatic differentiation to pull out the coeffs of the
        % operator:
        deriv = chebfun2deriv( op );
        
        % This is used to discretize the linear constrains:
        [bcrow, bcvalue] = constructbc(bcArg, bcpos,...
            een, bcn, dom, scl, order);
        
        % This is used to discretize the PDE:
        [CC, rhs, bb, gg, Px, Py, xsplit, ysplit] =...
            constructDiscretisation(N, f, m, n, flag);
        
        % Convert all the different user inputs for bc into uniform format:
        bc = createbc(bcArg, ends);
        
        % Method for deciding how to solve the matrix equation:
        X = denseSolve(N, f, m, n);
        
        % Remove trailing coefficients.
        a = truncate(a, tol);
        
    end
end

function N = constructor( N, varargin )
% Constructing the chebop2 object.

% Get chebfun2 preferences
prefs = chebfunpref();
tol = prefs.cheb2Prefs.eps;

% If empty input arguments then return an empty chebop2 object.
if ( isempty( varargin ) )
    return
end

% What domain is the operator defined on?
if ( numel( varargin ) > 1 )
    ends = varargin{2};      % second argument should be a domain.
    if ( length( ends ) == 4 )
        if ( (diff(ends(1:2)) > 0 )&& (diff(ends(3:4)) > 0 ) )  % valid domain?
            domain = ends;
        else
            error('CHEBOP2:CONSTRUCTOR:DOMAIN','Empty domain');
        end
    else
        error('CHEBOP2:CONSTRUCTOR:INPUT','Argument should be a domain given by four doubles')
    end
else
    if ( isa( varargin{1}, 'function_handle' ) )
        % pick the default domain
        rect1 = [-1,1];
        rect2 = [-1,1];
        domain = [rect1 rect2];
    elseif ( isa( varargin{1},'double' ) )
        N = chebop2(@(u) u, varargin{1});  % set up identity operator on the domain.
        return
    else
        error('CHEBOP2:INPUTS','First argument is not an operator or domain.')
    end
end


% First argument in the constructor is the operator.  If the operator is
% univariate then it's a constant coefficient PDE, otherwise assume it is a
% variable coefficient.

if ( isa(varargin{1},'function_handle') )
    op = varargin{1};
    
    if ( nargin(op) == 1 )     % The PDE has constant coefficients.
        
        %   OLD way! ( Use elimination )
        % Use an elimination procedure to extract out the constant
        % coefficients of the PDE.
        %         for jj = 0:maxorder
        %
        %             for kk = 0:maxorder
        %
        %                 % Test with the functions x^j*y^k because they are
        %                 % eliminated by terms containing j+1 and k+1 derivatives or
        %                 % higher.
        %                 const = factorial(jj) .* factorial(kk);
        %                 test = chebfun2(@(x,y) (x.^jj.*y.^kk)/const, domain);
        %
        %                 % By evaluating the result of op(f) at (0,0) we remove
        %                 % all the terms that contain j-1 and k-1 derivatives or
        %                 % less.
        %                 A(jj+1,kk+1) = feval(op(test),0,0);
        %
        %             end
        %
        %         end
        
        
        % Trust that the user has formed the chebfun2 objects outside of
        % chebop2.
        u = adchebfun2( chebfun2(@(x,y) x.*y, domain) );
        v = op( u );
        % If the PDO has constant coefficients then convert to double:
        try
            A = cell2mat(v.der.derCell).';
        catch
            % PDO has variable coefficients, keep them in a cell array:
            A = v.der.derCell;
        end
        
    elseif ( nargin(op) == 2 )
        error('Did you intend to have @(x,y,u)?')
    elseif ( nargin(op) == 3 )
        % The coefficients of the PDE are now variable coefficient.
        
        % Setup a chebfun2 on the right domain
        u = adchebfun2( chebfun2(@(x,y) x.*y, domain) );
        x = chebfun2( @(x,y) x, domain );
        y = chebfun2( @(x,y) y, domain );
        % apply it to the operator
        v = op( x, y, u );
        A = v.der.derCell;  % cell array of variable coefficients.
    else
        error('CHEBOP2:CONSTRUCTOR:INPUT',...
            'Operator should be @(u) or @(x,y,u).')
    end
else
    error('CHEBOP2:CONSTRUCTOR:INPUT',...
        'First argument should be an operator')
end

% Often the coefficients are obtained with small rounding errors
% and it is important to remove the very small non-zero ones to
% have rank(A) correct.
if iscell(A)
    for jj = size(A,1)
        for kk = size(A,2)
            if isa(A{jj,kk},'double') && abs(A{jj,kk})<10*tol
                A{jj,kk} = 0;
            end
        end
    end
else
    A( abs(A) < 10*tol ) = 0;
end

% Construct chebop2 object. The boundary conditions will be given later.
N.domain = domain;
N.op = op;
N.coeffs = A;

% Calculate xorder and yorder of PDE.
% Find the differential order of the PDE operator.
if ( iscell( A ) )
    xorder = size(A, 2) - 1;
    yorder = size(A, 1) - 1;
elseif ( min( size( A ) ) > 1 )
    xorder = find( sum( abs( A ), 2 ) > 100*tol, 1, 'last') - 1;
    yorder = find( sum( abs( A ) ) > 100*tol, 1, 'last' ) - 1;
else
    if ( size(A, 1) == 1 )
        yorder = length(A) - 1;
        xorder = 0;
    else
        xorder = length(A) - 1;
        yorder = 0;
    end
end
N.xorder = xorder;
N.yorder = yorder;

end