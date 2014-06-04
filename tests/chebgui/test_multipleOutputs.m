function pass = test_multipleOutputs(pref)
% This test ensures that the STRINGPARSER class is doing the correct thing when
% called with multiple outputs.

%% BVPS

% Strings we want to test.
str = { 'u(-1) = 1'
    'u(0)=1, u(2) + x = 3'
    'u = 1,v = 0,w = 3'
    'x + diff(sin(v))'
    'sum(u,0,.5) = 0'
    'fred(sin(x-y),u) = fred(cos(x-z),u)'};

% Store the number of strings
n = numel(str);
anFun = cell(n, 1);
indVarNames = anFun; varNames = anFun; pdeVarNames = anFun;
eigVarNames = anFun;
commaSeparated = zeros(n, 1);

for k = 1:n
    [anFun{k}, indVarNames{k}, varNames{k}, pdeVarNames{k}, eigVarNames{k}, ...
        commaSeparated(k)] = stringParser.str2anon(str{k}, 'bvp');
end

correctAnFun = {
    'feval(u,-1)-1'
    'feval(u,0)-1;feval(u,2)+x-3'
    'u-1;v;w-3'
    'x+diff(sin(v))'
    'sum(u,0,.5)'
    'fred(@(x,y)sin(x-y),u)-fred(@(x,z)cos(x-z),u)'};

correctIndVar = {
    {'', ''}
    {'x', ''}
    {'', ''}
    {'x', ''}
    {'', ''}
    {'x', ''}
    };

correctVarNam = {
    {'u'}
    {'u'}
    {'u'; 'v'; 'w'}
    {'v'}
    {'u'}
    {'u'}
    };

correctCommaSep = [0; 1; 1; 0; 0; 0];

passBVP = zeros(n, 1);
for k = 1:n
    passBVP(k) = ( strcmp(anFun{k}, correctAnFun{k}) && ...
        all( strcmp(indVarNames{k}, correctIndVar{k}) ) && ...
        all( strcmp(varNames{k}, correctVarNam{k}) ) && ...
        isempty(pdeVarNames{k}) && isempty(eigVarNames{k}) && ...
        ( commaSeparated(k) == correctCommaSep (k) ) );
end

%% EIGS

str = {'v''''-lambda*v = 0'
    'u''''+u'' = lam*(u + u'') + x*u'};
n2 = numel(str);
anFun = cell(n2, 1);
indVarNames = anFun; varNames = anFun; pdeVarNames = anFun;
eigVarNames = anFun;
commaSeparated = zeros(n2, 1);

for k = 1:n2
    [anFun{k}, indVarNames{k}, varNames{k}, pdeVarNames{k}, eigVarNames{k}, ...
        commaSeparated(k)] = stringParser.str2anon(str{k}, 'eig');
end

correctAnFun = {{'diff(v,2)'
    'v'}
    {'diff(u,2)+diff(u)-x.*u'
    'u+diff(u)'}};

correctIndVar = {
    {'', ''}
    {'x', ''}
    };

correctVarNam = {
    {'v'}
    {'u'}
    };

correctEigNam = {
    {'lambda'}
    {'lam'}
    };

correctCommaSep = [0; 0];
passEIG = zeros(n2, 1);
for k = 1:n2
    passEIG(k) = ( strcmp(anFun{k}{1}, correctAnFun{k}{1}) && ...
        strcmp(anFun{k}{2}, correctAnFun{k}{2}) && ...
        all( strcmp(indVarNames{k}, correctIndVar{k}) ) && ...
        all( strcmp(varNames{k}, correctVarNam{k}) ) && ...
        isempty(pdeVarNames{k}) && ...
        strcmp(eigVarNames{k}, correctEigNam{k}) && ...
        ( commaSeparated(k) == correctCommaSep (k) ) );
end

%% PDES

str = {'u_t+x*u'''' = u'''};
n3 = numel(str);

anFun = cell(n3, 1);
indVarNames = anFun; varNames = anFun; pdeVarNames = anFun;
eigVarNames = anFun;
commaSeparated = zeros(n2, 1);

for k = 1:n3
    [anFun{k}, indVarNames{k}, varNames{k}, pdeVarNames{k}, eigVarNames{k}, ...
        commaSeparated(k)] = stringParser.str2anon(str{k}, 'pde');
end
%%
correctAnFun = {'-x.*diff(u,2)+diff(u)'};
correctIndVar = {{'x', 't'}};
correctVarNam = {{'u'}};
correctPdeNam = {{'u_t'}};
correctCommaSep = 0;

passPDE = zeros(n3, 1);
for k = 1:n3
    passPDE(k) = ( strcmp(anFun{k}, correctAnFun{k}) && ...
        all( strcmp(indVarNames{k}, correctIndVar{k}) ) && ...
        all( strcmp(varNames{k}, correctVarNam{k}) ) && ...
        strcmp(pdeVarNames{k}, correctPdeNam{k}) && ...
        isempty(eigVarNames{k}) && ...
        ( commaSeparated(k) == correctCommaSep (k) ) );
end

%% Concatenate all types
pass = [passBVP; passEIG; passPDE];

end