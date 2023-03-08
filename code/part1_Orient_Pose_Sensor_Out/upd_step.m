function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
   
    z_t;
    Ct=eye(6,15); %Order of matrix is 6x15 as position and orientation are only considered in this kalman filter part
R=eye(6)*0.001;  %believe in our sensor measurement
Kt= (covarEst*Ct')*pinv((Ct*covarEst*Ct')+R); %kalman gain calculation formula

%Implementng Update step
uCurr= uEst+ (Kt*(z_t - (Ct*uEst)));
covar_curr= covarEst- (Kt*Ct*covarEst);     %formula for kalman filter is used as the system is already linearised in process model


return;
end

