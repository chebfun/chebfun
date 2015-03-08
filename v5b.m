% This function allows for a user defined
% epslevel that is enforced at every call
% of HAPPINESSCHECK in both the CHEBTECH
% and TRIGTECH classes.
%
% Jared L. Aurentz, March 2015

function v5b(epsilon)

   % global variable for eps
   global myeps

   % no input
   if nargin == 0
      myeps = eps;

   % user defined epsilon
   else if nargin == 1
      myeps = epsilon;

   % too many inputs
   else
      error('Maximum of one input for v5b!')
   end

end
