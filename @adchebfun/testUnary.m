function pass = testUnary(funcList, tol)
%TESTUNARY  Test that ADCHEBFUN is doing the correct thing for values,
%   derivatives and lineary information. Used for most unary operators (e.g.
%   sin, cos, exp, ...).

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  This method proceeds as follows:
%
%  1. Loops through the list of functions to be tested, passed by test methods 
%     in the test/adchebfun folder.
%  2. For each function handle:
%       a) Check that correct values for the function are computed, using the
%          valueTesting() method. 
%       b) Check that the correct derivative is computed, using the
%          taylorTesting() method.
%       c) Check that the correct linearity information was computed, from the
%          LIN variable returned by the valueTesting() method.
%  3. Return the pass vector with all results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default tolerance
if ( nargin < 2 )
    tol = 1e-2;
end

% Initialise vector with pass information
pass = zeros(3, numel(funcList));

% Do the tests.
for k = 1:numel(funcList)   
    
    % Call the valueTesting method, which also returns linearity information
    [err, lin] = adchebfun.valueTesting(funcList{k});

    % First, check that the computed function values match what we expect
    pass(1, k) = ( err == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(funcList{k});
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );
    
    % Check that we received the correct linearity information
    pass(3, k) = ( lin == 0 );
    
end

end
