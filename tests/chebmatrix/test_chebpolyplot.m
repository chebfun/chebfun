function pass = test_chebpolyplot(pref)
% This test just ensures that chebmatrix chebpolyplot() does not crash.

% Create an invisible figure to plot on:
hfig = figure('Visible', 'off');

%% Setup 

d = [-2 2];                   % function domain
x = chebfun(@(x) x, d);         

%% Create a couple of chebmatrices and try to call chebpolyplot on them

% Only CHEBFUNS
cm1 = [sin(x-.32); cos(100*(x+.3)); 3*tanh(50*sin(x-.2))];
cm2 = [sin(x-.32), cos(100*(x+.3)); ...
    3*tanh(50*sin(x-.2)) 4*sin(40*(x+.25))];

% CHEBFUNS and doubles
cm3 = [sin(x-.32); cos(100*(x+.3)); 3];
cm4 = [sin(x-.32), 100; ...
    3 4*sin((x+.25))];

pass(1) = doesNotCrash(@() chebpolyplot(cm1));
pass(2) = doesNotCrash(@() chebpolyplot(cm2));
pass(3) = doesNotCrash(@() chebpolyplot(cm3));
pass(4) = doesNotCrash(@() chebpolyplot(cm4));


%%

close(hfig);

end

function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch
    pass = false;
end

end
