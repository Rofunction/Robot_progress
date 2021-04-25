% C is vector with (nxn), and n is the number of 
%    Degrees of freedom of(DOF) the robotic Manipulators
function Cor=Coriolis(q,qdot)  %(n*n)
grav=[0,0,0];
n=length(q);
Cor=zeros(n,n);
Csq=zeros(n,n);
for i=1:n
    QD=zeros(1,n);
    QD(i)=1; % qdÎªµ¥Î»Õó
    tau = rnedyn(q, QD, zeros(1,n), grav);
    Csq (:,i) = Csq(:,i) + tau';
end

for i=1:n
    for j=i+1:n
        QD = zeros(1,n);
        QD(i) = 1;
        QD(j) = 1;
        tau = rnedyn(q, QD, zeros(1,n), grav);
        Cor(:,j)=Cor(:,j)+(tau' - Csq (:,j)-Csq (:,i))*qdot(i)/2;
        Cor(:,i)=Cor(:,i)+(tau' - Csq (:,j)-Csq (:,i))*qdot(j)/2;
    end
end
Cor=Cor + Csq * diag(qdot);
end