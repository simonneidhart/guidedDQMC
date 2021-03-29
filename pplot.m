load e1 %et_nowalker
load e2
load e3
load e4
%load w1 %walker pos
%load m1 %mu sigma

et_theory = 7.1147525409684089;
et_start = 8;

figure(1)
ll = length(e1);
x = 1:ll;
average_e1 = zeros(1,100);
average_e2 = zeros(1,100);
average_e3 = zeros(1,100);
average_e4 = zeros(1,100);
i1 = 1;
i2 = 100;
for i=1:100
    average_e1(i) = mean(e1(i1:i2,1));
    average_e2(i) = mean(e2(i1:i2,1));
    average_e3(i) = mean(e3(i1:i2,1));
    average_e4(i) = mean(e4(i1:i2,1));
    i1 = i1 + 100;
    i2 = i2 + 100;
end 
y = 1:100;
%semilogy(x,abs(e1(:,1) - et_theory),x,0.1./sqrt(x/10),'LineWidth',1.2)
semilogy(y,abs(average_e1 - et_theory),y,abs(average_e2 - et_theory),y,abs(average_e3 - et_theory),y,abs(average_e4 - et_theory),y,0.1./sqrt(y/10),'LineWidth',1.2)
xlabel('Iterations $\times$ 100','Interpreter','latex')
ylabel('$|\bar{E}(\tau)-E_0|$','Interpreter','latex')
legend('unguided','$c=0.1$','$c=0.01$','$c=0.001$','$\propto 1/\sqrt{t/10}$','Interpreter','latex')

figure(2)
plot(x,e1(:,2),x,e2(:,2),x,e3(:,2),x,e4(:,2),'LineWidth',1.0)
%plot(x(:),e1(:,1) - et_theory,x(:),e2(:,1) - et_theory,x(:),e3(:,1) - et_theory,x(:),e4(:,1) - et_theory,'LineWidth',1.0)
xlabel('Iterations','Interpreter','latex')
ylabel('Number of walkers','Interpreter','latex')
xlim([0 ll])
ylim([0.5*10^4 3.0*10^4])
legend('unguided','$c=0.1$','$c=0.01$','$c=0.001$','Interpreter','latex')

% l = length(positions(:,1));
% a = [1,1.5,2.0,0.5,3.0,1.0];
%     
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

figure(4)
plot(x(1000:end),e1(1000:end,1) - et_theory,x(1000:end),e2(1000:end,1) - et_theory,x(1000:end),e3(1000:end,1) - et_theory,x(1000:end),e4(1000:end,1) - et_theory,'LineWidth',1.0)
xlabel('Iterations','Interpreter','latex')
ylabel('$E(\tau)-E_0$','Interpreter','latex')
xlim([1000 ll])
legend('unguided','$c=0.1$','$c=0.01$','$c=0.001$','Interpreter','latex')


function y = gauss(x,ms)
    y = 1/sqrt(2*pi*ms(2)^2)*exp(-(x-ms(1)).^2/(2*ms(2)^2));
end 




