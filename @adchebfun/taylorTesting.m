function [order1, order2] = taylorTesting(f, hMax, numOut, plotting)
% VALUETESTING  Test that ADCHEBFUN is constructing the correct derivatives
%
% Here:
%   F is a function handle
%   numOut is an optional argument, used for methods with more than one outputs
%       (in particular, ellipj)

% TODO: Describe algorithm.

% Default NUMOUT is 1
if ( nargin < 3 )
    numOut = 1;
end

% By default, plotting is off.
if ( nargin < 4 )
    plotting = 0;
end

% Seed random generator to ensure same values.
seedRNG(6179);

% Initial u we evaluate f at
u = adchebfun(0.1*rand(8,1)+.5,[-1 1]);

% Chebfun used to create a perturbation
p = chebfun(0.01*rand(8,1)+.05,[-1 1]);

%% Taylor testing computations

% Vectors for storing information about convergence
nDiff1 = zeros(hMax, numOut);
nDiff2 = zeros(hMax, numOut);

% Factor we decrease by
fact = 5;

if ( numOut == 1)
    
    % Evaluate f at u
    fu = f(u);
    
    % Extract its Frechet derivative
    fuJac = fu.jacobian;
    
    for hCounter = 1:hMax
        % Compute a perturbation function
        pert = fact^-(hCounter+1)*p;
        
        % Evaluate f at the function obtained by adding the perturbation to u
        fuPert = f(u+pert);
        
        % Compute information for Taylor testing.
        
        % Expect this to be O(h)
        nDiff1(hCounter) = norm(fuPert - fu);
        
        % Use derivative information as well. Expect this to be O(h^2)
        nDiff2(hCounter) = norm(fuPert - fu - fuJac*pert);
    end
    
elseif ( numOut == 3)
    
    % Evaluate the function at u
    [fu, gu, hu] = f(u);
    
    % Extract all Frechet derivatives
    fuJac = fu.jacobian;
    guJac = gu.jacobian;
    huJac = hu.jacobian;
    
    for hCounter = 1:hMax
        % Compute a perturbation function
        pert = fact^-(hCounter+1)*p;
        
        % Evaluate f at the function obtained by adding the perturbation to u
        [fuPert, guPert, huPert] = f(u+pert);
        
        % Compute information for Taylor testing. We are actually testing the
        % correctness of three things at once, we hope for values around 1 in
        % the vectors order1 and order2 below.
        
        % Expect this to be O(h)
        nDiff1(hCounter, :) = [norm(fuPert - fu), norm(guPert - gu), ...
            norm(huPert - hu)];
        
        % Use derivative information as well. Expect this to be O(h^2)
        nDiff2(hCounter, :) = [norm(fuPert - fu - fuJac*pert), ...
            norm(guPert - gu - guJac*pert), norm(huPert - hu - huJac*pert)];
    end
else
    error('CHEBFUN:ADCHEBFUN:taylorTesting', ...
        'Unexpected number of output arguments of method being tested')
end
%% Compute the orders of convergence
% Expect this to have entries close to 1
order1 = bsxfun(@rdivide, diff(log(nDiff1)), diff(log(fact.^-(2:hMax+1))'));
% Expect this to have entries close to 2
order2 = bsxfun(@rdivide, diff(log(nDiff2)), diff(log(fact.^-(2:hMax+1))'));

% We expect order1 to have values around 1 and order2 to have values around 2.
% Thus, we can simply take the minimum values of each row, since incorrect
% derivatives will give values around 1 in order2 as well
order1 = min(order1, [], 2);
order2 = min(order2, [], 2);

%% Plotting

if ( plotting == 1 )
    loglog(fact.^-(1:hMax),nDiff1,'-*'), hold on
    % Expected rate of convergence
    loglog([1 5^-20],[1 5^-20],'k--')
    
    loglog(fact.^-(1:hMax),nDiff2,'r-*'), hold on
    % Expected rate of convergence
    loglog([1 5^-10],[1 5^-20],'k--')
    set(gca,'XDir','reverse'), shg
    grid on
    shg
    axis equal
end

end