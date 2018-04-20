
%%%%% LET'S CALL IT A "CONTROL" FILE %%%%%
clear all; close all; clc;

%Instantiate storage index for record keeping
macrorun=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZE ALGORITHM %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%User Defined SimOpt Parameters
Reps = 100; %total number of algorithmic replciations
T_setting = 250; %total number of simulations allotted per replication
n_0 = 10;   %number of initial design points
B = 1;     %number of replications per design point, set to 1 for deterministic problem
B_n0_setting=1;


%%%%%%%%%%%%%%%%%%%%%%  Problem to Solve  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%uncomment problem to solve, or input your own
% problem                 = testProblem1d();
% problem                 = testProblem1d_grlee();
 problem                 = testProblem2d_camel6();
% problem                 = testProblem_Griewank_kd();
% problem                 = testProblem_Rosen_kd();
% problem                 = testProblem_Zakharov_kd();


%Crowded EI level set threshold
alpha_lvl_set = 0.05; %i.e. EIs within 5% of maxEI

%%parameters for the TR algorithm, user defined
epsilon=1e-12; %for finite differencing
%for RC test and TR control
eta0=.25; 
eta1=.75;
delta=.5;
gamma=1.2;


%%%%%%%%%%%%%%%%   GRID  %%%%%%%%%%%%%%%%
%load('design_grid(100000x17)'); %use this grid for 17 dimensional problems
n0_grid = 2500; 
design_grid = Cdf(problem,n0_grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters for cross validation
MAXITERCROSS=1; %maximum number of cross validation iteration you want to perform
problem.varF = 0;
alpha = .05;
MAXIMUM = size(design_grid,1);

%Instantiate 'record' keeping process
x_storage = NaN(Reps,T_setting-n_0+1,problem.dim);

%%Begin Looping of Algorithm Replications

for L=12222:12222+Reps-1
    
%Begin Algorithm Replication 
macrorun = macrorun +1;

%%RNG Initialization
state = L;
    rng(state);
    
% Instantiate count of simulation effort expended
sim_count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation Optimization Begins %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Run "Simulations" %%%
    %Observe the locations associated to n_0, and store values    
reps=ones(n_0,1);
[M,mY] = RepSampSimple(problem,n_0,reps); %create initial design points x_0 and responses
x_0=M;
problem.x0 = x_0;
problem.y0 = problem.fun(problem.x0);
problem.xTrain = zeros(size(design_grid,1),size(design_grid,2)); %instantiate an xTrain
problem.yTrain = zeros(size(problem.xTrain,1),1);


%%Perfom the cross-validation procedure
   [B,T,problem.yTrain,problem.estNoiseVar,rep_cur,counter,problem.xTrain] = ...
       crossValProc_determin(MAXIMUM,x_0,n_0,B_n0_setting,alpha,T_setting,...
       problem.fHandle,problem.varF,MAXITERCROSS);
   
   sim_count = sim_count+n_0;
   
   %update all sampled points
   all_y = problem.yTrain;
   all_x = problem.xTrain;
   
   %%%%% Create First Record in x* %%%%%
   [y_min, index_min] = min(all_y(1:sim_count));
   x_storage(macrorun, sim_count-n_0+1,:) = all_x(index_min,:);
   
   find_max = 0; %needed if PI is utilized
   
while (sim_count < T_setting)
    
%%%%% Call Global OK Model %%%%%
 %Pass in x.train and y.train, Estimate global meta-model
    model = OK_model_kd_nugget(problem.xTrain(1:counter,:),...
    problem.yTrain(1:counter),0,2);   
     
%%%%% Call Expected Improvement EI function
 %pass in model from global OK, call EI function
 %discover next point to be centroid (x0)

    %%%%%% Using PI/EI here (look at name of fucntion called) %%%%%%%
   [EI_values, Varextrinsic, y_fit] = EIcalc_kd_pred(design_grid,...
       problem.xTrain(1:counter,:),model,problem.yTrain(1:counter)); 
   
    %%% Random Selection of next centroid based upon weighted PI's %%%
    %%%%%%% 'PRICE IS RIGHT STYLE', closest while being under %%%%%%%%%%
  %  Bid_mark = rand .* ones(size(PI_values,1),1);
  %  Price_line = PI_values - Bid_mark;
  %  closest_bid = 1;
  %  for i=1:size(Price_line,1)
  %      if Price_line(i)>0
  %         if Price_line(i)<closest_bid
  %             closest_bid = Price_line(i);
  %             x0=design_grid(i,:);
  %         end
  %      end
  %  end
    
    
 %If no bids are positive 
    %if closest_bid == 1 && find_max <= 2
    %    y_fit_max = -1e2;
    %    for i=1:size(y_fit)
    %        if y_fit(i)>y_fit_max
    %            y_fit_max = y_fit(i);
    %            x0 = design_grid(i,:);
    %        end
    %    end
    %    find_max = find_max + 1;
    %else
    %if closest_bid == 1
    %    x0 = rand().*(problem.ub-problem.lb)+problem.lb;
    %end
    %end
    
 %%------------------------------------------------------------%%   
 %%%%%%%%%%   CROWDING BASED CENTROID SELECTION (EI)    %%%%%%%%%
    
    EI_range = max(EI_values) - min(EI_values);
    relative_top_percent_threshold = EI_range*alpha_lvl_set;
 
    %%% Selection of next centroid based upon maxEI %%%
   [maxEI, Index] = max(EI_values);
   all_maxEI_pointer = 1;
   clear all_maxEI_locations;
   clear crowding_distance;
   %% Find all locations in top EI alpha level se %%
   for i=1:size(EI_values,1)
       if EI_values(i) >= (maxEI - relative_top_percent_threshold)
           all_maxEI_locations(all_maxEI_pointer,:)= design_grid(i,:);
           maxEI_value(all_maxEI_pointer,1) = EI_values(i);
           all_maxEI_pointer = all_maxEI_pointer+1;
       end
   end
   
   %Determine crowding distance at all locations in alpha level set
   for i=1:size(all_maxEI_locations,1)
       for j = 1:problem.dim
           DimWise_Crowd(j) = min(abs((ones(size(all_x,1),1).*all_maxEI_locations(i,j))-all_x(:,j)));
       end
       crowding_distance(i) = sum(DimWise_Crowd);
       clear DimWise_Crowd;
   end
   
   %Find the candidate with "max EI" that has the largest crowding distance
   %(i.e. is the furthest from all sampled points)
   [max_crowd, Index] = max(crowding_distance);
       x0 = all_maxEI_locations(Index,:);
   
   %%%% MAKING USE OF pure EI, no crowding %%%%
   %[maxEI, Index] = max(EI_values);
   %x0 = design_grid(Index,:);

 %%%% Execute Trust Region Search %%%%   
    %Create the static structure 'a' which contains constraint structure
    %for trust region optimization soltuion to be within the trust region.
     for i = 1:(2*problem.dim)
        if i <= problem.dim
            for j = 1:problem.dim
                if i == j
                    a(i,j) = 1;
                else
                    a(i,j) = 0;
                end
            end
        else
            for j = 1:problem.dim
                if (i-problem.dim) == j
                    a(i,j) = -1;
                else 
                    a(i,j) = 0;
                end
            end
        end
    end
    
    
    %%% Initialize TR to delta0 %%%
    TR_Bounds = [x0-problem.lb, problem.ub-x0, (problem.ub-problem.lb)./30];     %%% ATTN: HARD CODED TR_size0 value @.1 probably need a better criteion based on support range of inputs 
    TR_size0=min(TR_Bounds)*2;                
    TR_size=TR_size0;

    
    %%%%%% Gradient Initialization %%%%%% 
    %Creates a sampling plan of locations needed for finite differencing (stored as x)   
    for i = 1:(2*problem.dim)
        if i <= problem.dim
            for j = 1:problem.dim
                if i == j
                    x(i,j) = x0(1,j)+epsilon; 
                else
                    x(i,j) = x0(1,j);
                end
            end
        else
            for j = 1:problem.dim
                if (i-problem.dim) == j
                    x(i,j) = x0(1,j)-epsilon; 
                else
                    x(i,j) = x0(1,j);
                end
            end
        end
    end    
    
    %'sample' current centroid and record data
    f0=problem.fun(x0);
        sim_count=sim_count+1;
        all_y(sim_count)=f0;
        all_x(sim_count,:)=x0;
        %% Sample Taken => Record x* %%
            [y_min, index_min] = min(all_y(1:sim_count));
            x_storage(macrorun, sim_count-n_0+1,:) = all_x(index_min,:);
    
    %'sample' finite differencing plan and record data
    for i = 1:(2*problem.dim)
        sim_count = sim_count + 1;
        all_y(sim_count) = problem.fun(x(i,:));
        all_x(sim_count,:) = x(i,:);
        %% Sample Taken => Record x* %%
            [y_min, index_min] = min(all_y(1:sim_count));
            x_storage(macrorun, sim_count-n_0+1,:) = all_x(index_min,:);
    end
        
    %% These Derivest Suite functions handle any dimension function %%
    G1=jacobianest(problem.fHandle, x0);
    H = hessian(problem.fHandle, x0);


    %%--------------------------------------------------------%%
    %%%%%%% Begin Main Step of Trust Region Minimization %%%%%%%
    %%--------------------------------------------------------%%
    number_of_microloops=0;    
    
 while (norm(G1)>.05) && (TR_size>0.005)
    number_of_microloops=number_of_microloops+1;
    %construct quadratic model or linear model
        m = @(s) f0+s'*G1'+.5*s'*H*s;
    %construct contraints for step size to be within TR  
        b=(TR_size/2).*ones(2*problem.dim,1);
    %optimize step size through minimization of quadratic model
        sk=fmincon(@(s) quadratic_model(s,f0,G1',H),x0,a,b);
    
    %construct rho as the RC "test statistic"  
        fk = problem.fun(x0+sk);
        rho=(f0-fk)/(quadratic_model(zeros(1,problem.dim),f0,G1',H)-quadratic_model(sk,f0,G1',H));
            sim_count=sim_count+1;
            all_y(sim_count)=fk;
            all_x(sim_count,:)=(x0+sk);
            %% Sample Taken => Record x* %%
                [y_min, index_min] = min(all_y(1:sim_count));
                x_storage(macrorun, sim_count-n_0+1,:) = all_x(index_min,:);     
  
    %Compare trust region ratio with random number for probability of new centroid
    test = rand; 
    if test>(max(abs(sk))/(TR_size/2))
        TR_size=0;
    end
    
    %execute RC testing and TR control 
   if(rho<eta0)
        %reject candidate stepsize/centroid
        x0=x0;
        TR_size=TR_size*delta;
    else
        if(eta0<rho<eta1)
            %low pass of RC test
            x0=x0+sk;
            valid_bound = [2*(x0-problem.lb),2*(problem.ub-x0),TR_size];
            TR_size = min(valid_bound);
        else
            %high pass of RC test
            x0=x0+sk;
            valid_bound = [2*(x0-problem.lb),2*(problem.ub-x0), TR_size*gamma];
            TR_size = min(valid_bound);
        end

   end   
    
   
    %reinitialize main step to account for possibility that x0 moved 
    if (rho >= eta0) && (TR_size>0)
        %Creates a sampling plan of locations needed for finite differencing (stored as x)   
        for i = 1:(2*problem.dim)
            if i <= problem.dim
                for j = 1:problem.dim
                    if i == j
                        x(i,j) = x0(1,j)+epsilon; 
                    else
                        x(i,j) = x0(1,j);
                    end
                end
            else
                for j = 1:problem.dim
                    if (i-problem.dim) == j
                        x(i,j) = x0(1,j)-epsilon; 
                    else
                        x(i,j) = x0(1,j);
                    end
                end
            end
        end    

        %Given that we moved our centroid we know what the new value is at
        %that point from the RC test, and have already stored that data
        f0 = fk;

        %'sample' finite differencing plan and record data
        for i = 1:(2*problem.dim)
            sim_count = sim_count + 1;
            all_y(sim_count) = problem.fun(x(i,:));
            all_x(sim_count,:) = x(i,:);
            %% Sample Taken => Record x* %%
                [y_min, index_min] = min(all_y(1:sim_count));
                x_storage(macrorun, sim_count-n_0+1,:) = all_x(index_min,:);
        end   

        %reinitialize gradient for testing conditions at start of loop
        G1=jacobianest(problem.fHandle, x0);
        H =hessian(problem.fHandle, x0);
    end    
    
    
 end
    
   counter = counter + 1;
    
   %%%% Update the global model with the last centroid sampled %%%
   % % recall that problem.train trains the Gaussian process model % % 
           problem.xTrain(counter,:)=x0;
           problem.yTrain(counter)=f0;       
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TERMINATION SEQUENCE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlswrite3('your .xlsx file location here',x_storage,'Sheet1','A2',3);


   