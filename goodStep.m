%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% goodStep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% alpha = goodStep(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,...
%                      fPrime1,f2,fPrime2) finds 
% a global minimizer alpha within the bracket [brcktEndpntA,...
%                             brcktEndpntB] of the cubic polynomial 
% that interpolates f() and f'() at alpha1 and alpha2. Here ...
% f(alpha1) = f1, f'(alpha1) = fPrime1, 
% f(alpha2) = f2, f'(alpha2) = fPrime2.

function alpha = goodStep(brcktEndpntA,brcktEndpntB,...
                        alpha1,alpha2,f1,fPrime1,f2,fPrime2)


% Find interpolating Hermite polynomial in the z-space, 
% where z = alpha1 + (alpha2 - alpha1)*z
coeff = interpolatCubic(alpha1,alpha2,f1,fPrime1,f2,fPrime2);

% Convert bounds to the z-space
zlb = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
zub = (brcktEndpntB - alpha1)/(alpha2 - alpha1);

% Make sure zlb <= zub so that [zlb,zub] be an interval
if zlb > zub
  [zub,zlb] = deal(zlb,zub); % swap zlb and zub
end

% Minimize polynomial over interval [zlb,zub]
z = globalMin(zlb,zub,coeff); 
alpha = alpha1 + z*(alpha2 - alpha1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%