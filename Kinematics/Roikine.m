function thetaQ=Roikine(Poseref)
if length(Poseref)~=6
    error('input error, PoseRef is not 6 aixes')
end 
tic
ObtainQ=Roikine2(Poseref);
QSol=Roikine1(Poseref);
b=size(ObtainQ,1); 
c=size(QSol,1);
if b==c || b>c
    m=b;
elseif b~=c && c<b
    m=c;
end
for v=1:m
      poseError=Rofkine(ObtainQ(v,:)) - Poseref;
      poseError2=Rofkine(QSol(v,:)) - Poseref;
      if max(abs(poseError))>max(abs(poseError2)) && max(abs(poseError2))<10e-5
        thetaQ=QSol;
      elseif max(abs(poseError))<max(abs(poseError2)) && max(abs(poseError))<10e-5
          thetaQ=ObtainQ;
      end
end
toc
end 