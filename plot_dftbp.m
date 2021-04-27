load et_noWalkers
load walker_positions

figure(1)
ll = length(et_noWalkers);
x = 1:ll;
plot(x(1:end),(et_noWalkers(1:end,1)+5.6992283284),x(1:end),ones(1,ll)*0.0739272,'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('ZPE [Hartee]','Interpreter','latex')

zpe = mean(et_noWalkers(200:end,1))+5.6992283284

figure(2)
plot(x,et_noWalkers(:,2),'LineWidth',1.2)
xlabel('Iterations')
ylabel('Number of walkers')

nat = floor(length(walker_positions)/8);
s = 20.0;
c = [2,2,1,1,1,1,1,1];
c = repmat(c,1,nat);

figure(3)
scatter3(walker_positions(:,1),walker_positions(:,2),walker_positions(:,3),s,c,'filled')

err = zeros(800,1);
for i=1:800
    err(i) = mean(et_noWalkers(100:200+i,1))+5.6992283284 - 0.0739272;
end

figure(4)
plot(1:800,err,'LineWidth',1.2)
xlabel('Iterations')
ylabel('Error')
    
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





