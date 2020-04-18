function [prod] = CrossMe(v1, v2)
    % function to perform the cross product on symbolics
    prod = sym('prod', [3, 1]);
    prod(1) = v1(2)*v2(3) - v1(3)*v2(2);
    prod(2) = v1(3)*v2(1) - v1(1)*v2(3);
    prod(3) = v1(1)*v2(2) - v1(2)*v2(1);
end