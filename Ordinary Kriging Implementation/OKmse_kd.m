function y = OKmse_kd(X_0, X, sigma_z, Theta, Cov_Inverse)
% Calculate the MSE for the modified nugget effect model
% X_0 - all the possible design locations
% X - current design locations 
% sigma_z -the estimated process variance
% theta - the estimated sensitivity parameter
% cov_inverse - the inversed deterministic covariance matrix
% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

k = size(X,1);
temp = zeros(size(X_0,1),1);
Cov_Pred = zeros(k,1);
F = ones(k,1);  %%%F value

for c = 1:size(X_0,1)
    for i = 1:k
        %vecCheck = sum(X_0(c,:)-X(i,:));
        %vecCheckSquare=vecCheck.^2;
        Cov_Pred(i) = sigma_z.*exp(-(X_0(c,:)-X(i,:)).^2*Theta);
        %Cov_Pred(i) = sigma_z.*exp(-vecCheckSquare*Theta);
    end
    temp(c) = sigma_z.*( 1- Cov_Pred'*Cov_Inverse*Cov_Pred + ((1-F'*Cov_Inverse*Cov_Pred).^2)/(F'*Cov_Inverse*F));
    %temp(c) = sigma_z - Cov_Pred'*Cov_Inverse*Cov_Pred + ((1-F'*Cov_Inverse*Cov_Pred).^2)/(F'*Cov_Inverse*F);
end

y = temp;

end

