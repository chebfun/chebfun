% Test file for ADCHEBFUN trigonometric and related functions, part 4

function pass = test_trig4

% List of trigonometric functions to test.
trigFunctions = {@csch, @sec, @secd, @sech, @sin, @sind, @sinh, @tan, ...
    @tand, @tanh};

pass = zeros(2, numel(trigFunctions));

% Do the tests.
for k = 1:numel(trigFunctions)    
    % Call the taylorTesting method
    try
        % First, check that the computed function values match what we expect
        pass(1,k) = ( valueTesting(trigFunctions{k}) == 0 );
        
        [order1, order2] = taylorTesting(trigFunctions{k},4);
        % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be close
        % to 2.
        pass(2,k) = ( max(abs(order1-1)) < .01 ) & ( max(abs(order2-2)) < .01 );
    catch ME
        % Denote that an error happened when computing derivative, because the
        % version of development we branched off in Git didn't include the
        % required methods
        pass(1,k) = -1;
        pass(2,k) = -1;
    end
end

end