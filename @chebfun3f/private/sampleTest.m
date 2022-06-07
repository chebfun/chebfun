function happy = sampleTest(f, cf3Fhandle, tol, dom)


% Evaluate function
[x, y, z] = halton(30);
x = x*(dom(2)-dom(1))+dom(1);
y = y*(dom(4)-dom(3))+dom(3);
z = z*(dom(6)-dom(5))+dom(5);
for i = 1:30
    vFun(i) = f(x(i),y(i),z(i));
    vHandle(i) = cf3Fhandle(x(i),y(i),z(i));
end

% Perform tests
fail = (any(max(abs(vHandle - vFun)) > 10*tol));
happy = ~fail;
end


function [x, y, z] = halton(numpts)
p = haltonset(3,'Skip',1);
H = p(1:30,:);
x = H(:,1);
y = H(:,2);
z = H(:,3);
end