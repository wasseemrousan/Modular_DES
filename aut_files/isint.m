function [bool,idx] = isint(x)
    % Check whether input is integer or not
    % Inf and NaN are not integers
    if ~isnumeric(x)
       error('Input must be a numeric, not a %s.',class(x))
    end
    bool = (mod(x,1) == 0);
%     bool = (round(x) == floor(x));  % Other approach. Fails with Inf and
%     NaN.
    idx = find(bool); 
% Manolín Sept-2019
end