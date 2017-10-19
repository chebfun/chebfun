function pass = test_plotTree(~)
%TEST_PLOTTREE   Do TREEVAR computations, check that PLOTTREE works.

hfig = figure('Visible', 'off');

%% Basic computation:
u = treeVar();
v = cos(u);
w = sin(u);
t = v + w;

% Call both the static and non-static version of the plotting methods.
pass(1) = doesNotCrash(@() treeVar.plotTree(t.tree));
pass(2) = doesNotCrash(@() plot(t));

%% Introducing differentiation
u = treeVar();
myfun = @(u) 2 + diff(u,2);
s = myfun(u);
pass(3) = doesNotCrash(@() treeVar.plotTree(s.tree));
pass(4) = doesNotCrash(@() plot(s));

end


function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
