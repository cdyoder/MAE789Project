function [rotatedVect] = VectRotation(rotationMat,vect2Rotate)
    rotatedVect = simplify(rotationMat * vect2Rotate);
end