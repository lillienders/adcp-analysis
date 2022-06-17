function u = wake_function(eta,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    u = b(1)*log(eta)+ b(2) + b(3)*eta.^2.*(3-2*eta);
end

