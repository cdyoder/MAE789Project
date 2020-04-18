function [OvBO] = Transport1(OwB, rBO, t)
    % first derivative transport theorem
    OvBO = simplify(diff(rBO, t) + CrossMe(OwB, rBO));
end