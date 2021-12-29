%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFD.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradf = FFD(fun,xm,fm,hmin,hmax,typicalx)
%
% compute forward finite difference derivatives.

function gradf = FFD(fun,xm,fm,hmin,hmax,typicalx)

fm = fm(:);
n = numel(xm); 

xm = xm(:); typicalx = typicalx(:);
% Value of stepsize suggested in Trust Region Methods,
% Conn-Gould-Toint, section 8.4.3
CHG = sqrt(eps)*sign(xm).*max(abs(xm),abs(typicalx));
%
% Make sure step size lies within hmin and hmax
%
CHG = sign(CHG+eps).*min(max(abs(CHG),hmin),hmax);
variables = 1:n; gradf = zeros(n,1);
      
for gcnt=variables
   temp = xm(gcnt); xm(gcnt) = temp + CHG(gcnt);
   xOriginalShape = xm; xargin = xOriginalShape; fplus = fun(xargin);
   % Make sure it's in column form
   fplus = fplus(:); gradf(gcnt,1) = (fplus-fm)/CHG(gcnt); 
   xm(gcnt) = temp;
end % for 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





