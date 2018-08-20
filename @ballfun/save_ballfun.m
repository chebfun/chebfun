function save_ballfun(f,name)
% Save the tensor of coefficient of a ballfun in a mat file

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
save(name,"F");
end