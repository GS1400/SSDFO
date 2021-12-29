%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% globalMin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha = globalMin(lowerBound,upperBound,coeff) 
% finds a global minimizer alpha in the interval 
%                    lowerBound <= alpha <= upperBound of 
% the cubic polynomial defined by the coefficients in the 4-vector coeff. 

function alpha = globalMin(lowerBound,upperBound,coeff)

% Find stationary points 
stationaryPoint = roots([3*coeff(1) 2*coeff(2) coeff(3)]);

% Which among the two endpoints has a lower polynomial value?
[fmin,which] = min([polyval(coeff,lowerBound),polyval(coeff,upperBound)]);

if which == 1, alpha = lowerBound;else,alpha = upperBound;end

% If any of the stationary points is feasible, update the current
% global minimizer. If there's no stationary points, nothing is done
% below
if length(stationaryPoint) == 2 
  % Typical case: there are two stationary points. Either they are
  % both real or both are complex (with nonzero imaginary part) conjugate.
    if all(isreal(stationaryPoint))
        ok = lowerBound <= stationaryPoint(2)&& ...
              stationaryPoint(2) <= upperBound;
        if ok
            [fmin,which] = min([fmin,polyval(coeff,stationaryPoint(2))]);
            if which == 2, alpha = stationaryPoint(2);end
        end
        ok = lowerBound <= stationaryPoint(1) && ...
             stationaryPoint(1) <= upperBound;
        if ok
            [~,which] = min([fmin,polyval(coeff,stationaryPoint(1))]);
            if which == 2, alpha = stationaryPoint(1);end    
        end
    end
elseif length(stationaryPoint) == 1 % there is only one stationary point
    if isreal(stationaryPoint)   % the stationary point is real
        ok =  lowerBound <= stationaryPoint && ...
            stationaryPoint <= upperBound;
        if ok
            [~,which] = min([fmin,polyval(coeff,stationaryPoint)]);
            if which == 2, alpha = stationaryPoint; end
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
