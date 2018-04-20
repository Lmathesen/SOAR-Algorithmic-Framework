function [f mse]  = OK_predict(model,X_pred,regr_model,CovInverse)
% Build the modified nugget effect predictor based on the model given
% model - modified nugget effect kriging model, given by MNEK_model function 
% X_pred - locations to be predicted 
% regr_model - the underlying regression model for the mean function:
% regr_model = 0: constant mean function;
% regr_model = 1: linear mean function;
% regr_model = 2: quadratic mean function;

% Exmaple
%      M_predict  = MNEK_predict(model,X0,0);
% Using the parameter estimates of the MNEK model obtained from MNEK_model.m,
% the function predicts the response values at prediction points X0 with 
% a constant mean function

% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

% Obtain model parameters from MNEK_model
X = model.X;
min_X = model.min_X;
max_X = model.max_X;
[k d] = size(X);
theta = model.theta;
beta = model.beta;
beta_v = model.beta_v;
Z = model.Z;
Z_v = model.Z_v;
Z_m = model.Z_m;
L = model.L;
DZ_m = model.DZ_m;
D_L = model.D_L;
sigma_z = model.sigma_z;
corr_model = model.corr;

% get the size of the locations to be predicted
K = size(X_pred,1);
% get the regression model for the locations to be predicted
regr_pred = OK_regr(X_pred,regr_model);

% calculate distance matrix for prediction points
X_pred = (X_pred - repmat(min_X,K,1)) ./ repmat(max_X-min_X,K,1);


if corr_model == 2
    distXpred =  abs(repmat(reshape(X_pred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K])).^2;
else
    distXpred =  abs(repmat(reshape(X_pred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K]));
end

% calculate correlations between prediction points and design points
D = distXpred;
if corr_model == 3
    T = repmat(reshape(theta,[1 d 1]),[k 1 K]);
    R_pred = sigma_z*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
else
    R_pred = sigma_z*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
R_pred = reshape(R_pred,[k K 1]);

% calculate responses and MSE at prediction points 
f = regr_pred*beta + R_pred'*(L'\Z);

%mse = sigma_z.*ones(K,1) - diag(R_pred'*(D_L'\DZ_m)*R_pred) + sigma_z.*diag((R_pred'*(D_L'\DZ_m)*ones(size(R_pred))-eye(K))*(R_pred'*(D_L'\DZ_m)*ones(size(R_pred))-eye(K)))/(sum(sum(D_L'\DZ_m)));