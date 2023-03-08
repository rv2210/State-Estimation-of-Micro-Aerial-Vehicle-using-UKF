clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime,proj2Data] = init(datasetNum);

Z = sampledVicon(1:6,:);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.1*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0; %last time step in real time
pos = proj2Data.position;
pose = proj2Data.angle;
tic
for i = 1:length(sampledTime)
   
%%getting required data
acc= sampledData(i).acc;                    %Acceleration data extracted from student data to be given as input to prediction step function
dt= sampledTime(i)-prevTime;                %Sample time value extracted from student data to give as input to prediction step function previous time should be subtracted to get discrete time steps for the model
angVel=sampledData(i).omg;                  %Angular velocity data extracted from student data
zt= [transpose(pos(i,:)); transpose(pose(i,:))];    %z_t value taken from student data to be given as input to update step function
z_t= zt;

%% caling prediction and update step functions
[covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt);
[uCurr,covar_curr] = upd_step(z_t,covarEst,uEst);

savedStates(:,i)= uCurr; % current mean data is stored to savedstates to plot the data
uPrev= uCurr;               % making current mean to prev mean to be used for next mean calculation
covarPrev=covar_curr;       % making current covariance to prev covariance to be used for next covariance calculation
prevTime=sampledData(i).t;  % making current sample time as previous time to be used for next iteration 

end
toc
plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);