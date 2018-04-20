function f = cv_kd_OK(x, r, alpha)

    n = size(x,1);
    f = zeros(n,1);
 
    for i = 1:n
        x_1 = [x(1:i-1,:);x(i+1:n,:)];
        r_1 = [r(1:i-1);r(i+1:n)];
        model_1 = OK_model_kd_nugget(x_1,r_1,0,2);
        y_1 = OK_predict(model_1,x(i,:),0);
        cov_design = zeros(n-1);
        distance2=zeros(1,2);
        %compute the distance as the sum of squares on each dimension
        k=size(x_1,1);
        d=size(x_1,2);
        D_X = zeros(k, k, d);
        tempD_X = zeros(k*k,d);
        for h = 1:d
            hh=1;
            for ll=1:k
                for l = 1:k 
                    tempD_X(hh,h) = (x_1(ll,h) - x_1(l,h)).^2;
                    D_X(ll,l,h) = (x_1(ll,h) - x_1(l,h)).^2;
                    hh=hh+1;
                end
            end
        end
        distElem=zeros(1,2);
        hh=1;
        for j = 1:n-1
            for k = 1:n-1
                distElem=tempD_X(hh,:);
                cov_design(j,k) = model_1.sigma_z*exp(-distElem*model_1.theta);
                hh=hh+1;
            end
        end
        s_1 = OKmse_kd(x(i,:),x_1,model_1.sigma_z,model_1.theta,inv(cov_design));
        sigma = (r(i)-y_1)/s_1;
        if abs(sigma)>alpha
            f(i) = 1;
        end
    end
    
end