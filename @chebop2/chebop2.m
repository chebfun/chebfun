%CHEBOP2 constructor 
% 
% N = CHEBOP2(@(u) op(u)) constructs an operator N representing the
% operator given in @(u)op(u) acting on functions of two variables on
% [-1,1] by [-1,1]. 
%
% N = CHEBOP2(@(u) op(u),[a b c d]) constructs an operator N acting on
% functions of two variables defined on [a,b] by [c,d]. 
% 
% N = CHEBOP2(@(x,y,u) op(x,y,u),...) constructs a variable coefficient PDE
% operator. 

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


classdef chebop2
    
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
        %coeffcell={};   % cell array for storing variable coefficients
    end
    
    methods
        
        function N = chebop2(varargin)
            N = ctor(N, varargin{:});
        end
        
    end

    methods ( Static = true )
        X = BartelsStewart(A,B,C,D,E,xsplit,ysplit);
        bc = createbc(bcArg, ends);
        T = spconvermat(n,lam,k);
        M = MultMat(a,bn,varargin);
        diffMat = spdiffmat(n,k,varargin);
        [bcrow,bcvalue]=constructbc(bcArg,bcpos,een,bcn,dom,scl,order);
        deriv= chebfun2deriv( op )
    end
end

function N = ctor(N, varargin)
% Constructing the chebop2 object.

% Get chebfun2 preferences
prefs = chebfunpref(); 
tol = prefs.cheb2Prefs.eps;

% maximum differential order in the PDE.
% maxorder = 2;

% Set up empty arrays for coefficients. 
% C = cell(maxorder+1);
% A = zeros(maxorder+1); 

% If empty input arguments then return an empty chebop2 object.
if ( isempty(varargin) )
    return;
end

% What domain is the operator defined on?
if ( numel(varargin) > 1 )
    ends = varargin{2};      % second argument should be a domain.
    if ( length(ends) == 4 )
        if ( (diff(ends(1:2)) > 0 )&& (diff(ends(3:4)) > 0 ) )  % valid domain?
            domain = ends;
        else
            error('CHEBOP2:CONSTRUCTOR:DOMAIN','Empty domain');
        end
    else
        error('CHEBOP2:CONSTRUCTOR:INPUT','Argument should be a domain given by four doubles')
    end
else
    if isa(varargin{1},'function_handle') 
        % pick the default domain
        rect1 = [-1,1]; 
        rect2 = [-1,1];
        domain = [rect1 rect2];
    elseif isa(varargin{1},'double')
        N = chebop2(@(u) u, varargin{1});  % set up identity operator on the domain. 
        return;
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
    
%         %OLD way! % 
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
        
      try 
        %  NEW way! % 
        A = chebop2.chebfun2deriv(op); 
        A = rot90(A,2);
      catch 
          % Trust that the user has formed the chebfun2 objects outside of
          % chebop2. 
          u = adchebfun2(chebfun2(@(x,y) x.*y, domain)); 
          v = op(u); 
          A = v.der.derCell;
      end
    elseif (nargin(op)==2 )
        error('Did you intend to have @(x,y,u)?')
    elseif ( nargin(op) == 3 ) % The coefficients of the PDE are now variable coefficient.
        
        % setup a chebfun2 on the right domain 
        u = adchebfun2(chebfun2(@(x,y) x.*y, domain)); 
        x = chebfun2(@(x,y) x, domain); 
        y = chebfun2(@(x,y) y, domain);
        % apply it to the operator 
        v = op(x,y,u); 
        A = v.der.derCell;  % cell array of variable coefficients. 
        
        
%         for jj=0:maxorder
%             
%             for kk=0:maxorder
%                 % Test with the function x^j*y^k to eliminate terms with
%                 % derivatives of j+1, k+1 and higher. Now we have variable
%                 % coefficients so we are left with a variable coefficient
%                 % to approximate.
%                 
%                 const = factorial(jj).*factorial(kk);
%                 test = chebfun2(@(x,y) (x.^jj.*y.^kk)/const, domain);
%                 
%                 % First evaluate on a 5 by 5 grid to see if the variable
%                 % coefficient is actually a constant.
%                 x = chebpts(5); F = zeros(5);
%                 for j = 1:length(x)
%                     for k = 1:length(x)
%                         F(j,k) = feval(op(x(j),x(k),test),0,0);
%                     end
%                 end
%                 
%                 % If the variable coefficient is likely to be a constant
%                 % then it is faster to check for this case.
%                 if ( norm(F - mean(mean(F))) < 10 * tol )
%                     A(jj+1,kk+1) = mean(mean(F));  % constant is the mean
%                     C(jj+1,kk+1) = {[]};
%                 else
%                     % We have no choice but to call the constructor to work
%                     % out the variable coefficient.
%                     A(jj+1,kk+1) = 1; % Elimination procedure to get coeffs.
%                     C(jj+1,kk+1) = {chebfun2(@(x,y) feval(op(x,y,test),0,0),'vectorize')};
%                 end
%             end
%         end


    else
        error('CHEBOP2:CONSTRUCTOR:INPUT','Operator should be @(u) or @(x,y,u).')
    end
else
    error('CHEBOP2:CONSTRUCTOR:INPUT','First argument should be an operator')
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
%N.coeffcell = C;

% Calculate xorder and yorder of PDE. 
% Find the differential order of the PDE operator.
if iscell(A) 
    xorder = size(A,1)-1; 
    yorder = size(A,2)-1; 
elseif ( min(size(A)) > 1 )
    xorder = find( sum( abs( A ), 2 ) > 100*tol, 1, 'last') - 1;
    yorder = find( sum( abs( A ) ) > 100*tol, 1, 'last' ) - 1;
else
    if ( size(A,1) == 1 )
        yorder = length(A) - 1; xorder = 0;
    else
        xorder = length(A) - 1; yorder = 0;
    end
end
N.xorder = xorder; 
N.yorder = yorder; 

end
