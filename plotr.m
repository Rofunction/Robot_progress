t=linspace(0,10,10000);
[qd,dqd,ddqd]=DesirTrj_gd(t);

% figure (1)
% plot(this.t,this.qlist(:,1),'g',this1.t,this1.qlist(:,1),'b',this2.t,this2.qlist(:,1),'c',t.',qd(:,1),'r');
% legend('this(pdbig)','this1(pds)','this2(smc)','qd');
% xlabel('t');
% ylabel('qd');
% figure (2)
% plot(this.t,this.qlist(:,2),'g',this1.t,this1.qlist(:,2),'b',this2.t,this2.qlist(:,2),'c',t.',qd(:,2),'r');
% legend('this(pdbig)','this1(pds)','this2(smc)','qd');
% xlabel('t');
% ylabel('qd');
% figure (3)
% plot(this.t,this.qlist(:,3),'g',this1.t,this1.qlist(:,3),'b',this2.t,this2.qlist(:,3),'c',t.',qd(:,3),'r');
% legend('this(pdbig)','this1(pds)','this2(smc)','qd');
% xlabel('t');
% ylabel('qd');
% figure (4)
% plot(this.t,this.qlist(:,4),'g',this1.t,this1.qlist(:,4),'b',this2.t,this2.qlist(:,4),'c',t.',qd(:,4),'r');
% legend('this(pdbig)','this1(pds)','this2(smc)','qd');
% xlabel('t');
% ylabel('qd');
% figure (5)
% plot(this.t,this.qlist(:,5),'g',this1.t,this1.qlist(:,5),'b',this2.t,this2.qlist(:,5),'c',t.',qd(:,5),'r');
% legend('this(pdbig)','this1(pds)','this2(smc)','qd');
% xlabel('t');
% ylabel('qd');
% figure (6)
% plot(this.t,this.qlist(:,6),'g',this1.t,this1.qlist(:,6),'b',this2.t,this2.qlist(:,6),'c',t.',qd(:,1),'r');
% legend('this(pdbig)','this1(pds)','this2(smc)','qd');
% xlabel('t');
% ylabel('qd');

figure (1)
plot(this.t,this.qlist(:,1),'b',t.',qd(:,1),'r');
legend('actual','desired');
figure (2)
plot(this.t,this.qlist(:,2),'b',t.',qd(:,2),'r');
legend('actual','desired');
figure (3)
plot(this.t,this.qlist(:,3),'b',t.',qd(:,3),'r');
legend('actual','desired');
figure (4)
plot(this.t,this.qlist(:,4),'b',t.',qd(:,4),'r');
legend('actual','desired');
figure (5)
plot(this.t,this.qlist(:,5),'b',t.',qd(:,5),'r');
legend('actual','desired');
figure (6)
plot(this.t,this.qlist(:,6),'b',t.',qd(:,6),'r');
legend('actual','desired');



% figure (1)
% plot(this.t,this.qdlist(:,1),'b',t.',dqd(:,1),'r');
% figure (2)
% plot(this.t,this.qdlist(:,2),'b',t.',dqd(:,2),'r');
% figure (3)
% plot(this.t,this.qdlist(:,3),'b',t.',dqd(:,3),'r');
% figure (4)
% plot(this.t,this.qdlist(:,4),'b',t.',dqd(:,4),'r');
% figure (5)
% plot(this.t,this.qdlist(:,5),'b',t.',dqd(:,5),'r');
% figure (6)
% plot(this.t,this.qdlist(:,6),'b',t.',dqd(:,6),'r');

figure (7)
plot(this.t,this.e(:,1),'b');
legend('error 1');
figure (8)
plot(this.t,this.e(:,2),'b');
legend('error 2');
figure (9)
plot(this.t,this.e(:,3),'b');
legend('error 3');
figure (10)
plot(this.t,this.e(:,4),'b');
legend('error 4');
figure (11)
plot(this.t,this.e(:,5),'b');
legend('error 5');
figure (12)
plot(this.t,this.e(:,6),'b');
legend('error 6');



 














