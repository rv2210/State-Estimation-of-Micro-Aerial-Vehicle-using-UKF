function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
  
   g = [0; 0; -9.81];           %Defining gravity matrix and since gravity acts only in negative 
                                % z direction remaining elements are zero and this is in world frame
   roll = uEst(4);
   pitch = uEst(5);
   yaw = uEst(6);

   R=[cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch);
       cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll);
       -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];


   t_c_b = [-0.04; 0; -0.03];
   R_c_b = [0.7071 -0.7071 0; -0.7071 -0.7071 0 ; 0 0 -1]; % defining rotation matrix of camera wrt body from the images given in project 2
   R_b_c = transpose(R_c_b);
   T_c_b = [R_c_b, t_c_b; 0, 0, 0, 1];
   T_b_c = inv(T_c_b);

   %%defining parameters for UKF
    k= 1;
    alpha_val= 0.001;       %spread of sigma points about mean
    beta_val= 2;            % knowledge of gaussian, optimal value is 2 (reference textbook)
    n= 15;                  %dimensionality of x and q
    lambda= ((alpha_val^2)*(n+k)) - n;
    w0_m= lambda/(n+lambda);
    wi_m= 1/(2*(n+lambda));
    w0_c= (lambda/(n+lambda))+(1-(alpha_val^2)- beta_val);
    wi_c=1/(2*(n+lambda));

    Q_t =  0.03 * eye(3);
    %defining augmented matrix
    u_augt = uEst;
    covar_augt= covarEst;
    
    Z0_t = (R_c_b * transpose(R) * u_augt(7:9)) - (R_c_b * skew(T_b_c(1:3, 4)) * R_b_c * z_t(4:6));

    Z_t = Z0_t;
    sigma_pts = u_augt;


    cholesky_decomposition = chol(covar_augt,"lower"); %implementing cholesky decomposition

    %%computing sigma points    
    for i = 1 : n

        positiv_case = u_augt + (sqrt(n + lambda) * cholesky_decomposition(:,i));   %positive case means considered + in +/- case
        negativ_case = u_augt - (sqrt(n + lambda) * cholesky_decomposition(:,i));   %negative case means considered - in +/- case

        roll_pcase = positiv_case(4);
        pitch_pcase = positiv_case(5);
        yaw_pcase = positiv_case(6);
        roll_ncase = negativ_case(4);
        pitch_ncase = negativ_case(5);
        yaw_ncase = negativ_case(6);

        %calculating rotation matrices for sigma points
        R_positve=[cos(pitch_pcase)*cos(yaw_pcase), cos(yaw_pcase)*sin(pitch_pcase)*sin(roll_pcase) - cos(roll_pcase)*sin(yaw_pcase), sin(roll_pcase)*sin(yaw_pcase) + cos(roll_pcase)*cos(yaw_pcase)*sin(pitch_pcase);
           cos(pitch_pcase)*sin(yaw_pcase), cos(roll_pcase)*cos(yaw_pcase) + sin(pitch_pcase)*sin(roll_pcase)*sin(yaw_pcase), cos(roll_pcase)*sin(pitch_pcase)*sin(yaw_pcase) - cos(yaw_pcase)*sin(roll_pcase);
           -sin(pitch_pcase), cos(pitch_pcase)*sin(roll_pcase), cos(pitch_pcase)*cos(roll_pcase)];

        R_negative=[cos(pitch_ncase)*cos(yaw_ncase), cos(yaw_ncase)*sin(pitch_ncase)*sin(roll_ncase) - cos(roll_ncase)*sin(yaw_ncase), sin(roll_ncase)*sin(yaw_ncase) + cos(roll_ncase)*cos(yaw_ncase)*sin(pitch_ncase);
           cos(pitch_ncase)*sin(yaw_ncase), cos(roll_ncase)*cos(yaw_ncase) + sin(pitch_ncase)*sin(roll_ncase)*sin(yaw_ncase), cos(roll_ncase)*sin(pitch_ncase)*sin(yaw_ncase) - cos(yaw_ncase)*sin(roll_ncase);
           -sin(pitch_ncase), cos(pitch_ncase)*sin(roll_ncase), cos(pitch_ncase)*cos(roll_ncase)];

        Z_tpositive = (R_c_b * transpose(R_positve) * positiv_case(7:9)) - (R_c_b * skew(T_b_c(1:3,4)) * R_b_c * z_t(4:6));
        Z_tnegitive = (R_c_b * transpose(R_negative) * negativ_case(7:9)) - (R_b_c * skew(T_b_c(1:3,4)) * R_b_c * z_t(4:6));

        Z_t = [Z_t, Z_tpositive, Z_tnegitive];
        sigma_pts = [sigma_pts, positiv_case, negativ_case];

    end

    for i = 1 : (2*n + 1)

      if i == 1

          z = w0_m * Z_t(:,i);
      

       else
        z = z + (wi_m * Z_t(:,i));

       end

    end
%% calcultaing c_t and s_t uncertainity matrices

    for i = 1 : (2*n + 1)

        if i == 1

          c_t = w0_c * (sigma_pts(1:15 ,i) - uEst) * transpose(Z_t(:, i) - z);
          s_t = w0_c * (Z_t(:, i) - z) * transpose(Z_t(:, i) - z);

        else

          c_t = c_t + (wi_c * (sigma_pts(1:15 ,i) - uEst) * transpose(Z_t(:, i) - z));
          s_t = s_t + (wi_c * (Z_t(:, i) - z) * transpose(Z_t(:, i) - z));

        end

    end
%% calculating current mean and current covariance 
    s_t = s_t + Q_t;

    K_t = c_t/(s_t);  %calculating kalman gain

    uCurr = uEst + (K_t * (z_t(1:3) - z));
    covar_curr = covarEst - (K_t * s_t * transpose(K_t));






end

