function param = validateParameter(param, paramList, paramDesc)
    % Checks that the given parameter index or name is valid in the given list of parameters, then returns the
    % parameter name from the list.
    %
    % Examples
    % --------
    % >> param = validateParameter('d2o', {'D2O', 'SMW'});
    % This should return D2O as name check is case insensitive.
    %
    % Parameters
    % ----------
    % param
    %     The parameter name or index to validate. 
    % paramList
    %     The list of parameter to check for param.
    % paramDesc
    %     The description or type of the parameter to find in the list.
    arguments
        param
        paramList
        paramDesc {mustBeTextScalar} = 'Parameter'
    end

    paramList = cellstr(paramList);
    if isnumeric(param) && all(mod(param, 1) == 0)
        if (param < 1) || (param > length(paramList))
            throw(exceptions.indexOutOfRange(sprintf('%s %d is out of range (1 - %d).', paramDesc, param, length(paramList))));
        else
            param = paramList{param};
        end
    elseif isText(param)
        found = strcmpi(param, paramList);
        if ~any(found)
            throw(exceptions.nameNotRecognised(sprintf('Unrecognised %s name "%s".', lower(paramDesc), param)));
        end
        param = paramList{find(found, 1)};
    else
        throw(exceptions.invalidType(sprintf('Input %s should be a text or whole number.', lower(paramDesc))));
    end
end