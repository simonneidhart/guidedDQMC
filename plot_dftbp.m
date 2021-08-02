load et_noWalkers
load walker_positions

%e_min = -102.915322308400;
%e_min = -27.5893786053000;
e_min = -0.440620475522069;

figure(1)
ll = length(et_noWalkers);
x = 1:ll;
start = 50;
plot(x(start:end),(et_noWalkers(start:end,1)-e_min),'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('$E_T - E_{min}$ [Hartee]','Interpreter','latex')

zpe = mean(et_noWalkers(start:end,1))-e_min
st = std(et_noWalkers(start:end,1)-e_min)

figure(2)
plot(x,et_noWalkers(:,2),'LineWidth',1.2)
xlabel('Iterations')
ylabel('Number of walkers')

nat = 23;
s = 20.0;

wp = walker_positions;

c = ones(nat,1);
c(1) = 2;
c(2) = 2;
c = repmat(c,1,et_noWalkers(end,2));

%figure(3)
%scatter3(wp(:,1),wp(:,2),wp(:,3),s,'b','filled')

le = 1900;
start = 100;
err = zeros(le,1);
for i=1:le
    err(i) = mean(et_noWalkers(start:start+i,1))-e_min;
end

figure(4)
plot(1:le,err,'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('$Mean(E_T) - E_{min}$ [Hartee]','Interpreter','latex')
    
% x = -5:0.1:5;

% figure(3)
% subplot(3,2,1)
% hold on
% histogram(positions(:,1),50,'Normalization','pdf')
% plot(x,gauss(x,mu_sigma(1,:)),'LineWidth',1.5)
% xlabel('x_1')
% ylabel('walkers')
% hold off
% subplot(3,2,2)
% hold on
% histogram(positions(:,2),50,'Normalization','pdf')
% plot(x,gauss(x,mu_sigma(2,:)),'LineWidth',1.5)
% xlabel('x_2')
% ylabel('walkers')
% hold off
% subplot(3,2,3)
% hold on
% histogram(positions(:,3),50,'Normalization','pdf')
% plot(x,gauss(x,mu_sigma(3,:)),'LineWidth',1.5)
% xlabel('x_3')
% ylabel('walkers')
% hold off
% subplot(3,2,4)
% hold on
% histogram(positions(:,4),50,'Normalization','pdf')
% plot(x,gauss(x,mu_sigma(4,:)),'LineWidth',1.5)
% xlabel('x_4')
% ylabel('walkers')
% hold off
% subplot(3,2,5)
% hold on
% histogram(positions(:,5),50,'Normalization','pdf')
% plot(x,gauss(x,mu_sigma(5,:)),'LineWidth',1.5)
% xlabel('x_5')
% ylabel('walkers')
% hold off
% subplot(3,2,6)
% hold on
% histogram(positions(:,6),50,'Normalization','pdf')
% plot(x,gauss(x,mu_sigma(6,:)),'LineWidth',1.5)
% xlabel('x_6')
% ylabel('walkers')
% hold off
% 
% figure(4)
% ll = length(et_noWalkers);
% x = 1:ll;
% plot(x(100:end),abs(et_noWalkers(100:end,1) - et_theory),'LineWidth',1.2)
% xlabel('Iterations')
% ylabel('$|E(\tau)-E_0|$','Interpreter','latex')
% title('Error plot 2')
% legend('$|E(\tau)-E_0|$','Interpreter','latex')





