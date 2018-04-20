function [B,T,y,v,rep_cur,counter,x] = crossValProc_determin(MAXIMUM,x_0,n_0,B_n0,alpha,T,fhandle,noise_f,MAXITERCROSS)
%CROSSVALPROC Summary of this function goes here
%   Detailed explanation goes here
    iterCROSS=1;
    cross_OK=1;
    %Start adding point until cross validate the model or max number of cv
    %iteration's been performed
    while (cross_OK>0)
        x=x_0;
        y = zeros(MAXIMUM,1); % sample mean vector
        v = zeros(MAXIMUM,1); % sample variance vector
        rep_cur = zeros(MAXIMUM,1); % replications vector
        counter = n_0; % current number of design points

        %%% Obtaining observations for initial design
        for j = 1:counter
            y(j) = fhandle(x_0(j,:)); 
            rep_cur(j) = B_n0;
            %for k=1:B_n0
                %observations(j,k)=T_D(k,1);
            %end
        end

        cross_val = cv_kd_OK(x(1:counter,:),y(1:counter),alpha); % cross validation
        if (sum(cross_val)==0 || iterCROSS>=MAXITERCROSS)
            cross_OK=0;
            'cross-validation passed'
            B = B_n0;
        else
            iterCROSS=iterCROSS+1;
            B_n0=B_n0+10;
            B=B_n0;
            T=T+10*n_0;
        end   
    end

end

