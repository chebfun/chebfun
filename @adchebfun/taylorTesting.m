function [order1, order2, nDiff2] = taylorTesting(f, hMax, numOut, plotting)
%TAYLORTESTING   Test that ADCHEBFUN is constructing the correct derivatives.
%
% This method is called as follows:
%   [ORDER1, ORDER2, NDIFF2] = TAYLORTESTING(F, HMAX, NUMOUT, PLOTTING)
%
% The inputs are:
%   F        --  A function handle.
%   HMAX     --  (optional) Number of times we want to decrease the amplitude of
%                the perturbation used for testing. Default value: 4.
%   NUMOUT   --  (optional) Used for methods with more than one outputs (in
%                particular, ellipj). Default value: 1.
%   PLOTTING --  (optional) Value of 1 denotes that the results are to be
%                plotted (useful for debugging purposes). Default value: 0.
%
% The outputs are:
%   ORDER1   --  Computed values that correspond which should give first order
%                of convergence
%   ORDER2   --  Computed values that correspond which should give second order
%                of convergence
%   NDIFF2   --  Norm of difference between first order perturbation using AD
%                computed derivatives.
%
% See also: TAYLORTESTINGBINARY, VALUETESTING, VALUETESTINGBINARY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The following is taken from Birkisson's DPhil thesis [TODO: Citation?].
%
% Taylor testing is based on the observation that given an operator $F$, an
% argument $u$ and a perturbation to the argument $\delta u$, then by Taylor's
% theorem, one expects
%
%       T_1(h) := || F(u+h\delta u) - F(u) || -> 0 at O(h).
%
% However, if additionally, one has the \Frechet derivate of $F$ at $u$,
% $F'(u)$, then Taylor's theorem dictates that
%
%       T_2(h) := || F(u+h\delta u) - F(u) - h F'(u)\delta u || -> 0 at O(h^2).
%
% $T_1$ and $T_2$ are known as first- and second-order Taylor remainders,
% respectively.
%
% The procedure of Taylor testing is as follows. Let an operator $F$ be given.
% Assume that we have another operator $A$, and want to test whether $A$ is the
% \Frechet derivative of $F$. Then, to test whether the relation $A = F'$ holds,
% we repeatedly half $h$, and check whether
%
%       T_2(h/2)/T_2(h) ~ 1/4.

%% Parse inputs and initialise

if ( nargin < 2 )
    % Default hMax:
    hMax = 4;
end
if ( nargin < 3 )
    % Default NUMOUT is 1.
    numOut = 1;
end
if ( nargin < 4 )
    % By default, plotting is off.
    plotting = 0;
end

% Seed random generator to ensure the same values.
seedRNG(6179);

% Length of test function and perturbation:
N = 8;

% Initial u we evaluate f at:
u = adchebfun(chebfun(0.1*rand(N, 1) + .5, [-1 1]));

% Chebfun used to create a perturbation:
p = chebfun(0.01*rand(N, 1) + .05, [-1 1]);

%% Taylor testing computations

% Vectors for storing information about convergence:
nDiff1 = zeros(hMax, numOut);
nDiff2 = zeros(hMax, numOut);

% Factor we decrease by:
fact = 1/5;

% Evaluate F at u:
[fu{1:numOut}] = f(u);
for hCounter = 1:hMax
    % Compute a perturbation function:
    pert = fact^(hCounter+1)*p;

    % Evaluate f at the function obtained by adding the perturbation to u:
    [fuPert{1:numOut}] = f(u + pert);

    % Compute information for Taylor testing.

    % Expect this to be O(h)
    nDiff1(hCounter, :) = cellfun(@(u, v) norm(u - v), fuPert, fu);
    
    % Use derivative information as well. Expect this to be O(h^2)   
    nDiff2(hCounter, :) = ...
        cellfun(@(u, v) norm(u - v - v.jacobian*pert), fuPert, fu);
    
end


%% Compute the orders of convergence

% Expect this to have entries close to 1
order1 = bsxfun(@rdivide, diff(log(nDiff1)), diff(log(fact.^(2:hMax+1))'));
% Expect this to have entries close to 2
order2 = bsxfun(@rdivide, diff(log(nDiff2)), diff(log(fact.^(2:hMax+1))'));

% We expect order1 to have values around 1 and order2 to have values around 2.
% Thus, we can simply take the minimum values of each row, since incorrect
% derivatives will give values around 1 in order2 as well
order1 = min(order1, [], 2);
order2 = min(order2, [], 2);

%% Plotting
if ( plotting == 1 )
    loglog(fact.^(1:hMax), nDiff1, '-*'), hold on
    loglog(fact.^(1:hMax), nDiff2, 'r-*')
    % Expected rates of convergence
    loglog([1 fact^10], [1 fact^20], 'k--')
    loglog([1 fact^20], [1 fact^20], 'k--'); hold off
    set(gca, 'XDir', 'reverse'), grid on, axis('equal'), shg
end

end

