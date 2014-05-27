function f = diffx(f,k)

if nargin < 2 
    k = 1; 
end

f = diff(f,k,2);

end