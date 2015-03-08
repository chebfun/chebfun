% V5B  Chebfun v5 with prescribed epslevel.
%   This function forces the epslevel returned
%   by HAPPINESSCHECK in both the CHEBTECH
%   and TRIGTECH classes to be Matlab EPS.
function v5b

   % global variable for eps
   global myeps
   myeps = eps;
   
end
