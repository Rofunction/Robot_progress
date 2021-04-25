%{
        opts:   1) pose
                2) desired Motion
                3) noise // three parameters [tp,fat,peak], and time limit is [0,10];
                4) track
                5) tau
%}
function PlotFcn(model,opt)
    try
    qlist   = model.qlist;
    qdlist  = model.qdlist;
    endTime = model.endTime;
    t       = model.t;
    end 
    
    try
    qlist   = model.qlist;
    qdlist  = model.qdlist;
    endTime = model.endTime;
    t       = model.t;
    end 
    speed = floor(0.01*length(qlist));
% noise paratemers [tp, fat,peak]
    try
    tp  =model.tp;
    fat =model.fat;
    peak=model.peak;
    end
    
   ur5_mode; 
   
    if strcmp(opt,'pose')
        pose=zeros(length(qlist),6);
        for i=1:length(qlist)
           pose(i,:)=Fkine(qlist(i,:),'xyz'); 
        end
        [px,py,pz]=deal(pose(:,1),pose(:,2),pose(:,3));
        for i=1:speed:length(qlist)
            ur5.plot(qlist(i,:))        % 调用ur5去画轨迹（这里ur5只会显示运动，但不会显示实际轨迹）
            hold on;
            plot3(px(1:i),py(1:i),pz(1:i),'r-','linewidth',1);  %（显示实际轨迹）
            hold off;
        end
    end
    
    if strcmp(opt,'track')
       [px_d, py_d, pz_d] =DesiredTrj(t); 
       pose=zeros(length(qlist),6);
       for i=1:length(qlist)
           pose(i,:)=FK(qlist(i,:),'xyz');
       end
       [px,py,pz] =deal(pose(:,1),pose(:,2),pose(:,3));
       for i=1:speed:length(qlist)
       ur5.plot(qlist(i,:))
       hold on
       plot3(px_d(1:i),py_d(1:i),pz_d(1:i),'g-','linewidth',1.5);
       plot3(px(1:i)  ,py(1:i),  pz (1:i), 'r-','linewidth', 1 );
       hold off
       end
    end
    if strcmp(opt,'tau')
        tau=zeros(length(qlist),6);
       for i=1:length(qlist)
            M_q=(Inertia(qlist(i,:)));   
            Cor=Coriolis(qlist(i,:),qdlist(i,:));     
            G_q=(Grav(qlist(i,:),[0,0,-9.81])).';
            tau(i,:)=tau_t(t(i),qlist(i,:).',qdlist(i,:).',M_q,Cor,G_q).';
       end
       plot(t,tau)
    end
%     if strcmp(opt,'sim_results')
        
end

     
       
       
        
            
        
        
    
    
    
