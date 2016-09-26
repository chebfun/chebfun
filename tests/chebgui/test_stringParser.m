function pass = test_stringParser(pref)
% This test ensures that the STRINGPARSER class is doing the correct thing in a
% few different situations.

%% BVPS

% Strings we want to test.
str = { 'u(0) = 1'
        'u(-1) = 1'
        'u(0)=1, u(2) = 3'
        'u = 1,v = 0,w = 3'
        'u(end),u(1,left)'
        'diff(sin(u))'
        'diff(sin(u),2)'
        'diff(u^2) + sin(x)*u = 1'
        'feval(u,sqrt(2))'
        'feval(u,0,left)'
        'sum(u) = 0'
        'sum(u,0,.5) = 0'
        'volt(sin(x-y),u) = 0'
        'fred(sin(x-y),u) = fred(cos(x-z),u)'
        'u + fred(sin(2*pi*(y-x)),u) = feval(u,0,left)'};

% Store the number of strings
n = numel(str);
outBVP = cell(n, 1);

for k = 1:n
    outBVP{k} = stringParser.str2anon(str{k}, 'bvp');
end

correctBVP = {  '@(u) feval(u,0)-1'
                '@(u) feval(u,-1)-1'
                '@(u) feval(u,0)-1;feval(u,2)-3'
                '@(u,v,w) [u-1;v;w-3]'
                '@(u) feval(u,''end'');feval(u,1,''left'')'
                '@(u) diff(sin(u))'
                '@(u) diff(sin(u),2)'
                '@(x,u) diff(u.^2)+sin(x).*u-1'
                '@(u) feval(u,sqrt(2))'
                '@(u) feval(u,0,''left'')'
                '@(u) sum(u)'
                '@(u) sum(u,0,.5)'
                '@(x,u) volt(@(x,y)sin(x-y),u)'
                '@(x,u) fred(@(x,y)sin(x-y),u)-fred(@(x,z)cos(x-z),u)'
                '@(x,u) u+fred(@(x,y)sin(2.*pi.*(y-x)),u)-feval(u,0,''left'')'};

passBVP = zeros(n, 1);
for k = 1:n
    passBVP(k) = strcmp(outBVP{k}, correctBVP{k});
end

%% EIGS

str = {'u''''-lambda*u = 0'
       'u'' + sin(x)*u = lambda*u'
       'u''''+u'' = lambda*(u + u'') + u'};
n2 = numel(str);
outEIG = cell(n2, 1);

for k = 1:n2
    outEIG{k} = stringParser.str2anon(str{k}, 'eig');
end

correct = {{'@(u) diff(u,2)'
            '@(u) u'}
            {'@(x,u) diff(u)+sin(x).*u'
            '@(x,u) u'}
            {'@(u) diff(u,2)+diff(u)-u'
            '@(u) u+diff(u)'}};
        
passEIG = zeros(n2, 1);
for k = 1:n2
    for j = 1:numel(outEIG{k})
        passEIG(k) = strcmp(outEIG{k}{j}, correct{k}{j});
    end
end

%% PDES

str = {'u_t+x*u'''' = u'''};
n3 = numel(str);


outPDE = cell(n3, 1);
for k = 1:n3
    outPDE{k} = stringParser.str2anon(str{k}, 'pde');
end

correct = {'@(x,t,u) -x.*diff(u,2)+diff(u)'};

passPDE = zeros(n3, 1);
for k = 1:n3
    passPDE(k) = strcmp(outPDE{k}, correct{k});
end

%% Concatenate all types
pass = [passBVP; passEIG; passPDE];

end