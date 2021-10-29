load et_noWalkers
load walker_positions

% e_min = -102.915340231900; %h2@c60
%e_min = -27.5893786053000; %c12h10o
%e_min = -0.440620475522069; %peroskite
%e_min = -0.670556596200000; %h2
e_min = -102.2448840269; %c60
%e_min = -5.69922832840000; %c2h6
%unguided_res = 0.07386; 

figure(1)
ll = length(et_noWalkers);
x = 1:ll;
s = 100;
plot(x(s:end),(et_noWalkers(s:end,1)-e_min),'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('$E_T - E_{min}$ [Hartee]','Interpreter','latex')
title('guiding wf adjusted @ n=800')

zpe = mean(et_noWalkers(s:end,1))-e_min
st = std(et_noWalkers(s:end,1)-e_min)

figure(2)
plot(x,et_noWalkers(:,2),'LineWidth',1.2)
xlabel('Iterations')
ylabel('Number of walkers')

nat = 23;
s = 20.0;

wp = walker_positions;

c = ones(nat,3);
c(end,:) = 2;
c(end -1,:) = 2;
c = repmat(c,1,et_noWalkers(end,2));

figure(3)
scatter3(wp(:,1),wp(:,2),wp(:,3),s,'b','filled')
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])

le = 4500;
s = 500;
err = zeros(le,1);
for i=1:le
    err(i) = mean(et_noWalkers(s-99:s+i,1))- e_min;
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





