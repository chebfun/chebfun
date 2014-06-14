function pass = test_plotcoeffs(pref)
% This test just ensures that CHEBMATRIX/PLOTCOEFFS() does not crash.

% Create an invisible figure to plot on:
hfig = figure('Visible', 'off');

%% Setup 

d = [-2 2];                   % function domain
x = chebfun(@(x) x, d);         

%% Create a couple of CHEBMATRICES and try to call PLOTCOEFFS() on them

% Only CHEBFUNS
cm1 = [sin(x-.32); cos(100*(x+.3)); 3*tanh(50*sin(x-.2))];
cm2 = [sin(x-.32), cos(100*(x+.3)); ...
    3*tanh(50*sin(x-.2)) 4*sin(40*(x+.25))];

% CHEBFUNS and doubles
cm3 = [sin(x-.32); cos(100*(x+.3)); 3];
cm4 = [sin(x-.32), 100; ...
    3 4*sin((x+.25))];

pass(1) = doesNotCrash(@() plotcoeffs(cm1));
pass(2) = doesNotCrash(@() plotcoeffs(cm2));
pass(3) = doesNotCrash(@() plotcoeffs(cm3));
pass(4) = doesNotCrash(@() plotcoeffs(cm4));

%%

close(hfig);

end

function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME
    rethrow(ME)
    pass = false;
end

end
