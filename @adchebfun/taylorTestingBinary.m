function [order1, order2, nDiff2] = taylorTestingBinary(func, hMax, plotting)
%TAYLORTESTING   Test that ADCHEBFUN is constructing the correct derivatives
%   for binary operators.
%
% This method is called as follows:
%   [ORDER1, ORDER2, NDIFF2] = TAYLORTESTINGBINARY(F, HMAX, PLOTTING)
%
% The inputs are:
%   F        --  A function handle.
%   HMAX     --  (optional) Number of times we want to decrease the amplitude of
%                the perturbation used for testing. Default value: 4.
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
% See also: TAYLORTESTING, VALUETESTING, VALUETESTINGBINARY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Parse inputs and initialise

if ( nargin < 2 )
    % Default hMax:
    hMax = 3;
end
if ( nargin < 3 )
    % By default, plotting is off.
    plotting = 0;
end

% Seed random generator to ensure same values.
seedRNG(6179);

% Length of test functions and perturbations:
N = 8;

% Generate two arbitrary ADCHEBFUN objects to evaluate the function at:
u1 = adchebfun(chebfun(0.1*rand(N, 1) + .5, [-1 1]));
u2 = adchebfun(chebfun(0.1*rand(N, 1) + .5, [-1 1]));

% Create copies with re-seeded derivatives (i.e., 1x2 CHEBMATRIX derivatives)
v1 = seed(u1, 1, 2);
v2 = seed(u2, 2, 2);

% Also create two arbitrary CHEBFUNS to test with
w1 = chebfun(0.1*rand(N, 1) + .5);
w2 = chebfun(0.1*rand(N, 1) + .5);
 
% And two scalars to test with:
s1 = rand();
s2 = rand();
 
% Two CHEBFUNS used to create perturbations:
p1 = chebfun(0.01*rand(N, 1) + .05, [-1 1]);
p2 = chebfun(0.01*rand(N, 1) + .05, [-1 1]);

%% Create various combinations of variables

% Evaluate func at v1 and v2 -- ADCHEBFUN and ADCHEBFUN
v1v2 = func(v1, v2);

% Evaluate func at u1 and w2 -- ADCHEBFUN and CHEBFUN
u1w2 = func(u1, w2);

% Evaluate func at w1 and u2 -- CHEBFUN and ADCHEBFUN
w1u2 = func(w1, u2);

% Evaluate func at u1 and s2 -- ADCHEBFUN and SCALAR
u1s2 = func(u1, s2);

% Evaluate func at s1 and u2 -- SCALAR and ADCHEBFUN
s1u2 = func(s1, u2);

% Extract Frechet derivatives:
v1v2Jac = v1v2.jacobian;
u1w2Jac = u1w2.jacobian;
w1u2Jac = w1u2.jacobian;
u1s2Jac = u1s2.jacobian;
s1u2Jac = s1u2.jacobian;

%% Taylor testing computations

% Vectors for storing information about convergence
nDiff1 = zeros(hMax, 5);
nDiff2 = zeros(hMax, 5);

% Factor we decrease by
fact = 1/5;

for hCounter = 1:hMax
    % Compute a perturbation function
    pert1 = fact^(hCounter+1)*p1;
    pert2 = fact^(hCounter+1)*p2;

    % Evaluate the function handle passed in at the function obtained by adding
    % perturbations to the ADCHEBFUN objects
    v1v2Pert = func(v1 + pert1, v2 + pert2);
    u1w2Pert = func(u1 + pert1, w2);
    w1u2Pert = func(w1, u2 + pert2);
    u1s2Pert = func(u1 + pert1, s2);
    s1u2Pert = func(s1, u2 + pert2);
    
    % Compute information for Taylor testing.
    
    % Expect this to be O(h)
    nDiff1(hCounter,1) = norm(v1v2Pert - v1v2);
    nDiff1(hCounter,2) = norm(u1w2Pert - u1w2);
    nDiff1(hCounter,3) = norm(w1u2Pert - w1u2);
    nDiff1(hCounter,4) = norm(u1s2Pert - u1s2);
    nDiff1(hCounter,5) = norm(s1u2Pert - s1u2);
    
    % Use derivative information as well. Expect this to be O(h^2) for nonlinear
    % operators, and machine error for linear operators.
    
    % Need to store temporarily to be able to access the CHEBFUN block of the
    % resulting CHEBMATRIX.
    v1v2JacPert = v1v2Jac*[pert1; pert2];
    nDiff2(hCounter,1) = norm(v1v2Pert - v1v2 - v1v2JacPert{1});
    nDiff2(hCounter,2) = norm(u1w2Pert - u1w2 - u1w2Jac*(pert1));
    nDiff2(hCounter,3) = norm(w1u2Pert - w1u2 - w1u2Jac*(pert2));
    nDiff2(hCounter,4) = norm(u1s2Pert - u1s2 - u1s2Jac*(pert1));
    nDiff2(hCounter,5) = norm(s1u2Pert - s1u2 - s1u2Jac*(pert2));
end

%% Compute the orders of convergence

% Expect this to have entries close to 1
order1 = bsxfun(@rdivide, diff(log(nDiff1)), diff(log(fact.^(2:hMax+1))'));
% Expect this to have entries close to 2
order2 = bsxfun(@rdivide, diff(log(nDiff2)), diff(log(fact.^(2:hMax+1))'));

% We expect order1 to have values around 1 and order2 to have values around 2.
% Thus, we can simply take the minimum values of each row, since incorrect
% derivatives will give values around 1 in order2 as well
% order1 = min(order1, [], 2);
% order2 = min(order2, [], 2);

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
