clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.01*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0;
vel = proj2Data.linearVel;
angVel2 = proj2Data.angVel;
%% Calculate Kalmann Filter
tic
for i = 1:length(sampledTime)
    %% FILL IN THE FOR LOOP
     %% Fill in the FOR LOOP

%%getting required data
acc= sampledData(i).acc;                    %Acceleration data extracted from student data to be given as input to prediction step function
dt= sampledTime(i)-prevTime;                %Sample time value extracted from student data to give as input to prediction step function previous time should be subtracted to get discrete time steps for the model
angVel=sampledData(i).omg;                  %Angular velocity data extracted from student data
zt= [transpose(vel(i,:)); transpose(angVel2(i,:))];    %z_t value taken from student data to be given as input to update step function
z_t= zt;

%% caling prediction and update step functions
[covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt);
[uCurr,covar_curr] = upd_step(z_t,covarEst,uEst);

savedStates(:,i)= uCurr;    % current mean data is stored to savedstates to plot the data
uPrev= uCurr;               % making current mean to prev mean to be used for next mean calculation
covarPrev=covar_curr;       % making current covariance to prev covariance to be used for next covariance calculation
prevTime=sampledData(i).t;  % making current sample time as previous time to be used for next iteration 
end
toc
plotData(savedStates, sampledTime, sampledVicon, 2, datasetNum);