function [y] = ProblemF1d(x)
%PROBLEMF Summary of this function goes here
%   Detailed explanation goes here
%Specify the two-dimensional test function in this file
    temp = zeros(size(x,1),1);
    for i = 1:size(x,1)
       temp(i) = 10*(0.2*(x(i)-0.02)+1)*cos(13*(x(i)-0.02));
    end
    y = temp;
end

