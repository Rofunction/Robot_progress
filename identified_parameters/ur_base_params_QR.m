%{    
 -------------------------------------------------------------------------------
    Numerical Calculation of THE BASE INERTIAL PARAMETERS AND REGRESSOR MATRIX OF ROBOTS
    Note: 
        1¡¢ You should already have a function to compute the full regressor matrix of the
        robot. 
        2¡¢ This procedure has been verified.
        3¡¢ This procedure refers to M.CAUTIER's article 
            "Numerical Calculation of THE BASE INERTIAL PARAMETERS AND REGRESSOR
           MATRIX OF ROBOTS", 1990, IEEE.
--------------------------------------------------------------------------------
%}
% Seed the random number generator based on the current time
rng('shuffle');
% the limit parameters of robotic manipulator (ur5)  
q_min = -2*pi*ones(6,1);  
q_max = 2*pi*ones(6,1);
dq_max = deg2rad([180*ones(1,3),360*ones(1,3)]);
q2d_max = 2*ones(6,1); % it is chosen by us as it is not given in URDF

% ------------------------------------------------------------------------------
% Full standard dynamics parameters of the robot in symbolic form
% ------------------------------------------------------------------------------
m = sym ('m%d', [6,1], 'real');
h_x = sym ('h%d_x', [6,1], 'real');
h_y = sym ('h%d_y', [6,1], 'real');
h_z = sym ('h%d_z', [6,1], 'real');
ixx = sym ('i%d_xx',[6,1], 'real');
ixy = sym ('i%d_xy',[6,1], 'real');
ixz = sym ('i%d_xz',[6,1], 'real');
iyy = sym ('i%d_yy',[6,1], 'real');
iyz = sym ('i%d_yz',[6,1], 'real');
izz = sym ('i%d_zz',[6,1], 'real');

for i=1:6
    Id_qr_sym(:,i)=[ixx(i),ixy(i),ixz(i),iyy(i),iyz(i),izz(i),...
                       h_x(i),h_y(i),h_z(i),m(i)];
end
    [nLnkPrms, nLnks] = size (Id_qr_sym); 
    Id_qr_sym = reshape (Id_qr_sym , [nLnkPrms*nLnks ,1]);

% ------------------------------------------------------------------------------
% Find relation between independent columns and dependent columns
% ------------------------------------------------------------------------------
% Get observation matrix of full regressor matrix
Y_ob = [];
for i=1:20
    q_d = q_min + (q_max - q_min) .*rand(6,1);
    dq_d = -dq_max + 2*dq_max .*rand(6,1);
    ddq_d = -q2d_max + 2*q2d_max .*rand(6,1);  
    Y = full_regressor(q_d,dq_d,ddq_d);
    Y_ob = vertcat (Y_ob ,Y);
end
% QR decomposition with pivoting: Y_ob*E = Q*R
%   R is upper triangular matrix
%   Q is unitary matrix
%   E is permutation matrix
[Q,R,E]=qr(Y_ob);
rnk=rank(Y_ob);

% R = [R1 R2; 
%      0  0]
% R1 is bbxbb upper triangular and reguar matrix
% R2 is bbx(c-bb) matrix where c is number of standard parameters
R1=R(1:rnk,1:rnk);
R2=R(1:rnk,rnk+1:end);
beta=R1\R2;
beta(abs(beta)<sqrt(eps))=0; % get rid of numerical errors

% make sure that the relation Y_2=Y_1*beta holds
Y_1=Y_ob*E(:,1:rnk); 
Y_2=Y_ob*E(:,rnk+1:end);
if norm(Y_2-Y_1*beta)>1e-6
    fprintf('found the relationship between Y_1 with Y_2 is not correct\n');
    return;
end
% ------------------------------------------------------------------------------
% Find base parameters
% ------------------------------------------------------------------------------
X1=E(:,1:rnk)' * Id_qr_sym;     % independent variables
X2=E(:,rnk+1:end)' *Id_qr_sym;  % dependent variables
Id_qr_base = X1 + beta * X2;
%  or Id_qr_base = [eye(1:rnk,1:rnk) , beta]* E'*Id_qr_sym;
% Y_base = Y *E(:,1:rnk);

% ------------------------------------------------------------------------------
% Validation of obtained mappings
% ------------------------------------------------------------------------------
fprintf('Validation of mapping from standard parameters to base ones\n')
% Load full standard parameters
full_params=reshape(full_paramts,[nLnkPrms*nLnks, 1]);
 
% On random positions, velocities, aceeleations
 for i=1:100
    q_d = q_min + (q_max - q_min) .*rand(6,1);
    dq_d = -dq_max + 2*dq_max .*rand(6,1);
    ddq_d = -q2d_max + 2*q2d_max .*rand(6,1);
    
    Y_i = full_regressor(q_d,dq_d,ddq_d);
    tau_full = Y_i*full_params;
    Y_base = Y_i * E(:,1:rnk);
    Id_qr_base = [eye(rnk) , beta]* E'* full_params;
    tau_base = Y_base * Id_qr_base;
    
    num_errs(i)=norm(tau_full - tau_base);
 end
figure
plot (num_errs);
ylabel('||\tau - \tau_base||');
grid on
 
if ~all(num_errs<1e-6)
    fprintf('Validation failed');
    return;
end
    
% ------------------------------------------------------------------------------
% Create structure with the result of QR decompositon and save it
% for further use.
% ------------------------------------------------------------------------------   
base_QR = struct;
base_QR.numofbaseparameters = rnk;
base_QR.Id_qr_base = Id_qr_base;
base_QR.permutationMatrix = E; 

filename = 'base_QR.mat';
save (filename ,'base_QR');

    
    
    
    
    
    

