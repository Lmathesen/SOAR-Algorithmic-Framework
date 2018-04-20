function M_model = OK_model_kd_nugget(Xtrain, Ytrain, regr_model, corr_model) %, sigma_e)
% Build a Ordinary Kriging model for estimation
% X - design locations for the simulation inputs, size [k, d], k points with d
% dimensions 
% Y - observed simulation output values, size [k, 1], k points
% regr_model - the underlying regression model for the mean function:
% regr_model = 0: constant mean function;
% regr_model = 1: linear mean function;
% regr_model = 2: quadratic mean function;
% corr_model - the correlation model used for the spatial correlation
% corr_model = 0: linear correlation function
% corr_model = 1: exponential correlation function
% corr_model = 2: gaussian correlation function
% corr_model = 3: cubic spline correlation function

% Example
%       M_model = OK_model(X,Y, 0,1);
% This function uses OK model with a constant mean function and the 
% exponential correlation function to fit the data, (X,Y)
% where X is the design points and Y are the outputs at design points X 

% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

[k d] = size(Xtrain);
%k is the number of records in Xtrain

if (size(Ytrain,1) ~= k)
    error('Simulation outputs and inputs should have the same number of rows.')
elseif(size(Ytrain,2) ~= 1)
    error('Simulation outputs should only have one value at each location.')
end

%if (size(sigma_e,1) ~= k)
 %   error('Simulation outputs and inputs should have the same number of rows.')
%elseif(size(sigma_e,2) ~= 1)
 %   error('Simulation outputs should only have one value at each location.')
%end

if ((regr_model ~= 0)&(regr_model ~= 1)&(regr_model ~= 2))
    error('please type of regression model: 0 = constant, 1 = linear, 2 = quadratic.')
end

if ((corr_model ~= 0)&(corr_model ~= 1)&(corr_model ~= 2)&(corr_model ~= 3))
    error('please type of correlation function: 0 = linear, 1 = exponential, 2 = Gaussian, 3 = cubic spline.')
end

% RANDOM WHITE NOISE % sigma_e = diag(sigma_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    sigma_e = diag(sigma_e);
%data normalization
min_X = min(Xtrain);  
max_X = max(Xtrain);
%GG = repmat(min_X,k,1);
Xtrain = (Xtrain - repmat(min_X,k,1)) ./ repmat(max_X-min_X,k,1);

% calculate the distance between design locations
tmp = 0;
D_X = zeros(k, k, d);
tempD_X = zeros(k*k,d);
for h = 1:d
    hh=1;
    for i = 1:k
        for l = 1:k 
            
            tempD_X(hh,h) = (Xtrain(i,h) - Xtrain(l,h));
            
            D_X(i,l,h) = (Xtrain(i,h) - Xtrain(l,h));
            hh=hh+1;
        end
    end
end

%call the regression model
regr = OK_regr(Xtrain,regr_model);

% initialize parameters, least square estimation for beta and sample
% variance for sigma_z, 
beta_0 = (regr'*regr)\(regr'*Ytrain); %basic mean function when regr = 1 (constant)
sigma_z0 = var(Ytrain-regr*beta_0);   %variance of residuals (res = observ-mean)

%delta_0 = .5*sigma_z0;
%if (delta_0 > 1)
%    delta_0 = .99999;

theta_0 = zeros(d,1);                 %initialize hyperparameters, 1 per dimension
if (corr_model == 0 || corr_model ==3)  
    theta_0(1:d) = 0.5;               %if using linear or cubic correlation, initalize all theta to 0.5
else
                                      %o.w. theta is: log(2)/d * (mean absolute value of all distances)^(-2 or -1) 
    theta_0(1:d) = (log(2)/d)*(mean(abs(tempD_X)).^(-corr_model));
end

%initialize correlation matrix
% calculate the correlation matrix R based on the initial theta selected and correlation model selected
R = OK_corr(corr_model,theta_0,D_X);           %R=correlation matrix\

% initialize the lower bounds 
a=8;
delta_lb = max((max(eig(R))*(cond(R)-exp(a)))/(cond(R)*(exp(a)-1)),0);
%lob_delta =  max(delta_lb);
lob_sigma_z = 0.001*sigma_z0;    
lob_theta = 0.001*ones(d,1); 

lob = [lob_sigma_z;lob_theta];
% maximize profile log-likelihood function (-"logPL")
% subject to lower bounds on sigma_z and theta
myopt = optimset('Display','off','MaxFunEvals',1000000,'MaxIter',500);

%estimate hyperparameters
%estimate hyperparameters
try
    parms = fmincon(@(x) OK_lh_kd_nugget(x,k,d,D_X,Ytrain,regr,corr_model,delta_lb),[sigma_z0;theta_0],[],[],[],[],lob,[],[],myopt); 
catch
    fail = 1;
    count = 0;
    while fail == 1&& count < 300
        try
            parms = fmincon(@(x) OK_lh_kd_nugget(x,k,d,D_X,Ytrain,regr,corr_model,delta_lb),[rand();rand(length(theta_0),1)],[],[],[],[],lob,[],[],myopt); 
            fail = 0;
        catch
            fail = 1;
            count = count + 1;
        end
    end
end
% record MLEs for sigma_z and theta 
sigma_z = parms(1);
theta = parms(2:length(parms));


% calculate the correlation matrix R based on the theta and correlation model selected
R = OK_corr(corr_model,theta,D_X);           %R=correlation matrix\


CR = sigma_z*(R+delta_lb.*eye(size(R,1),size(R,2)));
DR = sigma_z*R;
[U0,pd0] = chol(CR);
%%[U1,pd1] = chol(sigma_z*R);
%%[U2,pd2] = chol(sigma_e); 
if(pd0>0)

    save data;
 %   error('covariance matrix is nearly singular');
    %end
end
CR=U0;
DR=U0;
%cholesky decomposition
L = U0';
D_L = U0';
L_inv = inv(L);
R_inv = L_inv'*L_inv;
beta = inv(regr'*R_inv*regr)*(regr'*(R_inv*Ytrain));
beta_v = inv(regr'*R_inv*regr)*(regr'*(R_inv));

% output MLEs and other things useful in prediction
M_model.sigma_z =  sigma_z;
M_model.min_X = min_X;
M_model.max_X = max_X;
M_model.regr =  regr;
M_model.beta = beta;
M_model.beta_v = beta_v;
M_model.theta = theta;
M_model.X = Xtrain;
M_model.corr = corr_model;
M_model.L = L;
M_model.D_L = D_L;
M_model.Z = L\(Ytrain-regr*beta);
M_model.Z_v = L\(eye(max(size(Ytrain)))-regr*beta_v);
M_model.Z_m = inv(L);
M_model.DZ_m = inv(D_L);