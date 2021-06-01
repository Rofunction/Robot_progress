%{
    Note: the condition number is 


%}
function out=traj_cond(opt_vars,traj_pars,base_QR)

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
     
% Compute trajectory (Fouruer series + fifth order polynomail)     
     [q,dq,ddq]=Fourier_series_trj(t,wf,A,B,q_init,N,C);
     
% the base permutation matrix is obtained by using QR decomposition
     E1=base_QR.permutationMatrix(:,1:base_QR.numofbaseparameters);
% Get the base regressor martix and observer martix
     W=[];
     for i=1:length(t)
         Y=full_regressor(q(:,i),dq(:,i),ddq(:,i))*E1;
         W=vertcat(W,Y);
     end
% computes the observer martix's condition number 
     out=cond(W);
end