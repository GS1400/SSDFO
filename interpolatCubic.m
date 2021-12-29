%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% interpolatCubic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coeff = interpolatCubic(alpha1,alpha2,f1,fPrime1,f2,fPrime2) 
% determines the coefficients of the cubic polynomial that 
% interpolates f and f' at alpha1 and alpha2; that is, 
% c(alpha1) = f1, c'(alpha1) = fPrime1, c(alpha2) = f2, 
% c'(alpha2) = fPrime2.


function coeff = interpolatCubic(alpha1,alpha2,f1,fPrime1,f2,fPrime2)


deltaAlpha = alpha2 - alpha1;
coeff(4) = f1;
coeff(3) = deltaAlpha*fPrime1;
coeff(2) = 3*(f2 - f1) - (2*fPrime1 + fPrime2)*deltaAlpha;
coeff(1) = (fPrime1 + fPrime2)*deltaAlpha - 2*(f2 - f1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%