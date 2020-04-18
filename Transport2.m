function [OaBO] = Transport2(OwB, OaB, rBO, t)
    % second derivative transport theorem
    OaBO = simplify(diff(rBO, t, t) + 2*CrossMe(OwB, diff(rBO, t)) + ...
        CrossMe(OaB, rBO) + CrossMe(OwB, CrossMe(OwB, rBO)));
end