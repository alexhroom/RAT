function chi2 = chiSquared(shiftedData,reflectivity,nParams)
% Chi-squared function is used to evaluate the goodness of fit.
% It is a measure of the difference between the observed and expected
% reflectivity.

    y = shiftedData(:,2);
    eb = shiftedData(:,3);
    fit = reflectivity(:,2);
    
    N = max(size(y));
    if N <= nParams
        N = (nParams + 1);
    end
    terms = ((y-fit)./eb).^2;
    n = find(terms == Inf);
    if ~isempty(n)
        terms(n) = 0;
    end
    chi2 = sum(terms)/(N-nParams);

end