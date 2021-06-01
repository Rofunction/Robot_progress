function [c,ceq]=traj_cnt(opt_vars,traj_pars,base_QR)

% Trajectory parameters
     N  = traj_pars.N;
     t  = traj_pars.t;
     T  = traj_pars.T;
     wf = traj_pars.wf;
     q_init = traj_pars.q_init;
% -------------------------------------------------------------------------------- 
% As paramters of the trajectory are in a signle vector we reshape them as
% to feed the function that computes the trajectory  
% --------------------------------------------------------------------------------
     AB=reshape(opt_vars,[12,N]);
     A=AB(1:6,:);  % sin coefficients
     B=AB(7:12,:); % cos coefficients
% --------------------------------------------------------------------------------
% Solve the coefficient of fifth order polynomial, the fifth order polynimial
% is added to fourier series in order to guarantee that positions, 
% velocities and accelerations are zero in the beginning and at time T.
% --------------------------------------------------------------------------------
     C=Solve_C(T,A,B,N,wf,q_init);
     
% codition number of the observer martix
% --------------------------------------------------------------------------------
% The function computes constraints on trajectory for trajectoty
% optimization needed for dynamic parameter identification more 
% precisely for searching feasibly solutions with condition number
% less than 100 that guarantees good estimate
% --------------------------------------------------------------------------------
     numcondition=traj_cond(opt_vars,traj_pars,base_QR);
     
% Compute trajectory (Fouruer series + fifth order polynomail)     
     [q,dq,ddq]=Fourier_series_trj(t,wf,A,B,q_init,N,C);
     
% Inequality constraint
     c(1:6)   = traj_pars.q_min - min(q,[],2);
     c(7:12)  = max(q,[],2) - traj_pars.q_max;
     c(13:18) = max(dq,[],2) - traj_pars.dq_max;
     c(19:24) = max(ddq,[],2) - traj_pars.ddq_max;
     
% more precisely for searching feasibly solutions
     c(25)= numcondition - 100;
     
% equality constraint
     ceq=[];
end