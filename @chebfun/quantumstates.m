function varargout = quantumstates(varargin)
%QUANTUMSTATES    Compute and plot Schroedinger eigenstates.
%   This program computes and plots eigenvalues lambda and eigenfunctions u
%   for the equation Lu = lambda*u, where L is the Schroedinger operator
%   defined by Lu(x) = -h^2*u"(x) + V(x)*u(x).  Here h is a small parameter
%   and the potential function V is given as a Chebfun. The domain of the
%   problem is the domain of V, with boundary conditions u=0 at both ends.
%
%   Inputs:
%
%       QUANTUMSTATES(V) plots 10 eigenstates for h=0.1
%       QUANTUMSTATES(V, n) plots n eigenstates for h=0.1
%       QUANTUMSTATES(V, h), h noninteger, plots 10 eigenstates for given h
%       QUANTUMSTATES(V, n, h) plots n eigenstates for given h
%       QUANTUMSTATES(..., 'noplot') produces no plot
%
%   Outputs:
% 
%       D = QUANTUMSTATES(...) returns a vector D of eigenvalues
%       [U, D] = QUANTUMSTATES(...) returns a quasimatrix U of eigenfunctions
%       and a diagonal matrix of eigenvalues
%
%   Examples:
%
%       x = chebfun('x', [-3, 3]);
%       V = x.^2;                 % harmonic oscillator, or
%       V = abs(x);               % absolute value, or
%       V = (x.^2-1).^4;          % double well

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
% Nick Trefethen, January 2012

%% Parsing of inputs:
noplot = strcmpi(varargin{nargin}, 'noplot');      % check for no plot
nargin1 = nargin;                                  % no of input args
if ( noplot )
    nargin1 = nargin1 - 1; 
end                    
V = varargin{1};                                   % potential function
n = 10;                                            % default no of states
h = 0.1;                                           % default constant
if ( nargin1 == 3 )
    n = varargin{2}; 
    h = varargin{3};
end
if ( nargin1 == 2 )
    v2 = varargin{2};
    if ( v2 == round(v2) )
        n = v2;
    else
        h = v2; 
    end
end

%% Eigenvalue computation:
[xmin, xmax] = domain(V);                          % domain of problem
% Create a CHEBOP with Dirichlet BCs
L = chebop([xmin, xmax]);
L.lbc = 0; L.rbc = 0;
L.op = @(x,u) -h^2*diff(u,2) + V.*u;               % Schroedinger operator
[U, D] = eigs(L, n, 'sr');                         % compute evals/efuns
d = diag(D);                                       % vector of evals
[d, ii] = sort(d);                                 % sort them
U = U(:,ii);                    

%% Outputs:
if ( nargout == 2 )
    varargout = {U, diag(d)};
else
    varargout = {d};
end

if ( noplot )
    % If we're not plotting, then we return.
    return
end

%% Plot:

holdState = ishold;

% Plot the potential function:
LW = 'linewidth';
plot(V, 'k', LW, 2, 'jumpline', '-k'), hold on   
s = sprintf('h = %4g      %d eigenstates', h, n);
title(s, 'fontsize', 12)

% Vertical limits:
ymax = max(d); 
ymin = min(V); 
ydiff = ymax - ymin; 
ymax = ymax + .2*ydiff; 
ymin = ymin - 0*ydiff;       

% V values at endpoints:
Vxmin = feval(V, xmin); 
Vxmax = feval(V, xmax);
if ( ymax > Vxmin )                              % The potential 
    plot(xmin*[1 1], [ymax, Vxmin], 'k', LW, 2)  %   V(x) effectively
end                                              %   goes to infinity
if ( ymax > Vxmax )                              %   at the endpoints,
    plot(xmax*[1 1], [ymax Vxmax], 'k', LW, 2)   %   so we make the
end                                              %   plot show this.
dx = .05*(xmax - xmin); 
dy = .25*ydiff/max(5, n);

% Plot the eigenfunction, lifted by the corresponding eigenvalue:
W = dy*U;
W = num2cell(W);
for j = 1:n
    umm = minandmax(W{j});
    if ( umm(2) < -umm(1) )
        W{j} = -W{j}; 
    end
    W{j} = W{j} + d(j);
end
W = horzcat(W{:});
plot(W, LW, 1.5)

% Plot V(x) again (so that black ends up on top):
plot(V, 'k', LW, 2, 'jumpline', '-k')
if ( ymax > Vxmin )
    plot(xmin*[1, 1], [ymax, Vxmin], 'k', LW, 2)
end                                          
if ( ymax > Vxmax )
    plot(xmax*[1, 1], [ymax Vxmax], 'k', LW, 2)    
end

% Set axis:
axis([xmin - dx, xmax + dx, ymin - dy, ymax]), drawnow

if ( ~holdState )
    % Stop holding:
    hold off
end
  
end                                            
