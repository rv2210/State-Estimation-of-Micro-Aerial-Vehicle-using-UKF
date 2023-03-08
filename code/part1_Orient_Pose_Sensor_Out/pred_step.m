function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 

%%calculating X_dot matrix
x= uPrev(1);
y= uPrev(2);
z= uPrev(3);
roll= uPrev(4);
pitch= uPrev(5);
yaw= uPrev(6);
vx= uPrev(7);
vy= uPrev(8);
vz= uPrev(9);
bg1= uPrev(10);
bg2= uPrev(11);
bg3= uPrev(12);
ba1= uPrev(13);
ba2= uPrev(14);
ba3= uPrev(15);
wx= angVel(1);
wy= angVel(2);
wz= angVel(3);
ax= acc(1);
ay= acc(2);
az= acc(3);
g= [0 ;0 ;-9.81]; %Defining gravity matrix and since gravity acts only in negative z direction remaining elements are zero and this is in world frame
  
bg=[bg1;bg2;bg3];   %Defining gyroscopic bias in x, y, z axis as bg1 bg2 bg3
ba=[ba1;ba2;ba3];   %Defining accelerometer bias in x, y, z axis as ba1 ba2 ba3
am=[ax;ay;az];      %defining accelerometer matrix
wm=[wx;wy;wz];      % defining angular velocity matrix
t=wm-bg;
u=am-ba;

%rotation matrix of euler angle parametrization Z-Y-X

R=[cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch); 
    cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll);
    -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];

G = [cos(pitch)*cos(yaw), -sin(yaw), 0;
    cos(pitch)*sin(yaw), cos(yaw), 0; 
    -sin(pitch), 0, 1] ;              

%defining elements of X_dot matrix

X1_dot=[vx;vy;vz];
X2_dot= pinv(G)*R*(t);
X3_dot= (g + (R * (u)));

X_dot= [X1_dot; X2_dot; X3_dot;0;0;0;0;0;0];


%Defining parameters for UKF prediction calculation
k= 1;
alpha_val= 0.001;   %spread of sigma points about mean
beta_val= 2;        % knowledge of gaussian, optimal value is 2 (reference textbook)
n_prime= 27;        %dimensionality of x and q, n'=n+ng=15+12
lambda_prime= ((alpha_val^2)*(n_prime+k)) - n_prime;
w0_m= lambda_prime/(n_prime+lambda_prime);%weight for mean calculation
wi_m= 1/(2*(n_prime+lambda_prime));
w0_c= (lambda_prime/(n_prime+lambda_prime))+(1-(alpha_val^2)+ beta_val); %weight for covariance calculation
wi_c=1/(2*(n_prime+lambda_prime));

%%defining augmented matrices 
X= X_dot;
Q_t=eye(12);

P_aug= [covarPrev zeros(15,12);
        zeros(12,15) Q_t];
u_aug=[uPrev; zeros(12,1)];
sigma_pt= u_aug;
cholesky_decomp= chol(P_aug,"lower");  %implementing cholesky decomposition

%%computing sigma points
for i= 1:n_prime
     
    positive_case= u_aug + ((sqrt(n_prime+lambda_prime)*cholesky_decomp(:,i))); %positive case means considered + in +/- case
    negative_case= u_aug - ((sqrt(n_prime+lambda_prime)*cholesky_decomp(:,i))); %negative case means considered - in +/- case
    
    %getting orienations of sigma points
    roll_pc= positive_case(4);
    pitch_pc= positive_case(5);
    yaw_pc= positive_case(6);
    roll_nc= negative_case(4);
    pitch_nc= negative_case(5);
    yaw_nc= negative_case(6);

    %computing rotation matrices of sigma points
       R_pc= [cos(pitch_pc)*cos(yaw_pc), cos(yaw_pc)*sin(pitch_pc)*sin(roll_pc) - cos(roll_pc)*sin(yaw_pc), sin(roll_pc)*sin(yaw_pc) + cos(roll_pc)*cos(yaw_pc)*sin(pitch_pc); 
                 cos(pitch_pc)*sin(yaw_pc), cos(roll_pc)*cos(yaw_pc) + sin(pitch_pc)*sin(roll_pc)*sin(yaw_pc), cos(roll_pc)*sin(pitch_pc)*sin(yaw_pc) - cos(yaw_pc)*sin(roll_pc);
                -sin(pitch_pc), cos(pitch_pc)*sin(roll_pc), cos(pitch_pc)*cos(roll_pc)];

       R_nc= [cos(pitch_nc)*cos(yaw_nc), cos(yaw_nc)*sin(pitch_nc)*sin(roll_nc) - cos(roll_nc)*sin(yaw_nc), sin(roll_nc)*sin(yaw_nc) + cos(roll_nc)*cos(yaw_nc)*sin(pitch_nc); 
            cos(pitch_nc)*sin(yaw_nc), cos(roll_nc)*cos(yaw_nc) + sin(pitch_nc)*sin(roll_nc)*sin(yaw_nc), cos(roll_nc)*sin(pitch_nc)*sin(yaw_nc) - cos(yaw_nc)*sin(roll_nc);
            -sin(pitch_nc), cos(pitch_nc)*sin(roll_nc), cos(pitch_nc)*cos(roll_nc)];
        
       
       G_pc = [cos(pitch_pc)*cos(yaw_pc), -sin(yaw_pc), 0;
                cos(pitch_pc)*sin(yaw_pc), cos(yaw_pc), 0; 
                 -sin(pitch_pc), 0, 1] ;
    
       G_nc = [cos(pitch_nc)*cos(yaw_nc), -sin(yaw_nc), 0;
                    cos(pitch_nc)*sin(yaw_nc), cos(yaw_nc), 0; 
                     -sin(pitch_nc), 0, 1] ;
    
       %% calculating Xt matrix and propagating sigma points through nonlinear function g
       Xt_pc = [positive_case(7:9);
                 (pinv(G_pc) * R_pc * (wm - positive_case(10:12) - positive_case(16:18))); 
                 (g + (R_pc * (am - positive_case(13:15) - positive_case(19:21))));
                  positive_case(22:27)];
    
       Xt_nc = [negative_case(7:9);
                   (pinv(G_nc) * R_nc * (wm - negative_case(10:12) - negative_case(16:18))); 
                   (g + (R_nc * (am - negative_case(13:15) - negative_case(19:21))));
                   negative_case(22:27)];
    
       
        
    
        X=[X, Xt_pc, Xt_nc];
        sigma_pt= [sigma_pt, positive_case, negative_case ];

end

    for i= 1:((2*n_prime)+1)

        Xt(:,i)= sigma_pt(1:15,i) + (dt*X(1:15,i));
    end


%%calculating predicted mean and predicted covariance 

   for i= 1:((2*n_prime)+1)
       if i==1
           uEst= w0_m*Xt(:,i);
           
       else
           uEst= uEst+(wi_m* Xt(:,i));
           
       end
   end
    for i= 1:((2*n_prime)+1)
       if i==1
          
           covarEst = (w0_c * (Xt(:,i) - uEst) * transpose(Xt(:, i) - uEst));
       else
           
           covarEst = covarEst + (wi_c * (Xt(:,i) - uEst) * transpose(Xt(:,i) - uEst));
       end
   end
   return

end



