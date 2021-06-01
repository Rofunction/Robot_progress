%{


%}
    if isempty(gcp('nocreate'))
        parpool;
    end
tic
    load('base_QR.mat');
% Choose optimization algorithm: 'patternsearch', 'GA'
    optalgorithm = 'patternsearch';

% Trajectory parameters 
    traj_pars.T      = 25;
    traj_pars.wf     = 2*pi/traj_pars.T;
    traj_pars.t_smp  = 2e-1;
    traj_pars.t      = 0:traj_pars.t_smp:traj_pars.T;
    traj_pars.N      = 7;            % number of harmonics
    traj_pars.q_init = deg2rad([0,-90,0,-90,0,0]');
    
% the limit parameters of robotic arm (ur5)
    q_min  = -2*pi*ones(6,1);
    q_max  =  2*pi*ones(6,1);
    dq_max = deg2rad( [180*ones(1,3),360*ones(1,3)]' );

% Trajectory limit for positions for safety
    traj_pars.q_min  = -deg2rad([180,180,100,180,90,90]');
    traj_pars.q_max  =  deg2rad([180,0,100,0,90,90]');
    traj_pars.dq_max = dq_max;
    traj_pars.ddq_max= [2,1,1,1,1,2.5]';

% ----------------------------------------------------------------------------------
%   Optimited: the purpose of optimization is to get the most suitable
%   parameters of A,B,C.
% ----------------------------------------------------------------------------------
% the Setting of patternsearch method
A=[]; b=[];Aeq=[];beq=[];lb = []; ub = [];

    if strcmp (optalgorithm,'patternsearch')
        x0 = rand(6*2*traj_pars.N,1);
%     x0 = reshape([A B], [6*2*traj_par.N, 1]);
        option_parh.optimoptions = 'patternsearch';
        option_parh.Display = 'iter';
        option_parh.ConstraintTolerance = 1e-6;
        option_parh.MaxFunctionEvaluations = 1e+6;
        option_parh.MaxTime = Inf;
        option_parh.FunctionTolerance = 10;
        option_parh.StepTolerance = 1e-1;
        
        [x,favl] = patternsearch (@(x)traj_cond(x, traj_pars, base_QR),...
                                     x0, A, b, Aeq, beq, lb, ub,...
                            @(x)traj_cnt(x, traj_pars, base_QR),option_parh); 
    elseif strcmp(optmznAlgorithm, 'ga')
        x0 = rand(6*2*traj_pars.N,1);
        option_ga = optimoptions('ga');
        option_ga.Display = 'iter';
        option_ga.PlotFcn = 'gaplotbestf'; % {'gaplotbestf', 'gaplotscores'}
        option_ga.MaxGenerations = 50;   
        option_ga.PopulationSize = 1e+3; % in each generation.
        option_ga.InitialPopulationRange = [-100; 100];
        option_ga.SelectionFcn = 'selectionroulette';
        
        [x,favl] = ga(@(x)traj_cond(x, traj_pars, base_QR), x0, A, b, Aeq,...
                            beq,lb,ub,@(x)traj_cnt(x, traj_pars, base_QR),...
                         option_ga);
    end
    
% ----------------------------------------------------------------------------------
%   Plot the optimized trajectory
% ----------------------------------------------------------------------------------  
    AB = reshape(x,[12,traj_pars.N]);
    A = AB(1:6,:); B = AB(7:12,:);
    C= Solve_C( traj_pars.T, A, B, traj_pars.N, traj_pars.wf, traj_pars.q_init );
    [q,dq,ddq]=Fourier_series_trj( traj_pars.t, traj_pars.wf, A,...
                                  B, zeros(6,1), traj_pars.N, C);
    
    figure
    subplot(3,1,1)
    plot(traj_pars.t,q);
    ylabel('$q$','interpreter','latex');
    grid on
    legend('q1','q2','q3','q4','q5','q6');
    subplot(3,1,2)
    plot(traj_pars.t,dq);
    ylabel('$\dot{q}$','interpreter','latex');
    grid on
    legend('dq1','dq2','dq3','dq4','dq5','dq6');
    subplot(3,1,3)
    plot(traj_pars.t,ddq);
    ylabel('$\ddot{q}$','interpreter','latex');
    grid on
    legend('ddq1','ddq2','ddq3','ddq4','ddq5','ddq6');

% ----------------------------------------------------------------------------------
%   save the parameters of optimized trajectory (A,B,C)
% ----------------------------------------------------------------------------------  
    pathToFolder = 'optmt_parameters';
    st = strcat('N',num2str(traj_pars.N),'T',num2str(traj_pars.T));
    if strcmp (optalgorithm,'patternsearch')
    filename = strcat(pathToFolder,'patternSearch_',st,'QR.mat');
    elseif strcmp(optalgorithm,'ga')
    filename = strcat(pathToFolder,'ga_',st,'.mat'); 
    end
    save(filename,'A','B','C','traj_pars');
toc


