function [pi_0, s_0, y_0] = PIcalc_kd_pred(x_0,x,M_model,y) % EI calculator

[curr_best, curr_best_ind] =  min(y);

cov_design = zeros(size(x,1));
k=size(x,1);
d=size(x,2);
D_X = zeros(k, k, d);
tempD_X = zeros(k*k,d);
for h = 1:d
    hh=1;
    for ll=1:k
        for l = 1:k 
            tempD_X(hh,h) = (x(ll,h) - x(l,h)).^2;
            D_X(ll,l,h) = (x(ll,h) - x(l,h)).^2;
            hh=hh+1;
        end
    end
end

distElem=zeros(1,2);
hh=1;
for j = 1:size(x,1)
    for k = 1:size(x,1)
        distElem=tempD_X(hh,:);
        cov_design(j,k) = M_model.sigma_z*exp(-distElem*M_model.theta)+0.0001;
        hh=hh+1;
    end
end

b_0 = ones(size(x_0,1),1); %design for constant mean regression?(vector of 1's as b)
y_0 = OK_predict(M_model,x_0,0,inv(cov_design)); 
s_0 = OKmse_kd(x_0,x,M_model.sigma_z,M_model.theta,inv(cov_design));

%%%% PLOT THE PREDICTED SURFACE %%%%
%[sortedX, sortIndex] = sort(x_0); %sort x_0 in ascending order
%sortedY = y_0(sortIndex); %sort y_0 maintaing original x_0/y_0 index pairs
%figure();
%plot(sortedX, sortedY) %plot of the global prediction over design_grid (x_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
found=0;
while (i<=size(x_0,1) && found==0)
    if x_0(i,:) == x(curr_best_ind,:)
        curr_best = y_0(i);
        found=1;
    else
        i=i+1;
    end
end
counts = size(x_0,1);
pi_0 = zeros(size(x_0,1),1);
s_0 = sqrt(s_0);
epsilon = 1e-2;

for i = 1:counts
    if s_0(i)>epsilon
        pi_0(i) =  normcdf((curr_best-y_0(i))/s_0(i),0,1);
    end
    if s_0(i)<=epsilon
        pi_0(i)=0;
    end
    
end

%%normalize probabilities
min_p = min(pi_0);
max_p = max(pi_0);
pi_0 = (pi_0 - repmat(min_p,size(x_0,1),1)) ./ repmat(max_p-min_p, size(x_0,1),1);

end
