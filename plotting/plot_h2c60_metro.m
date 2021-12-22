load et_noWalkers00
load et_noWalkers01
load et_noWalkers02
load et_noWalkers03
load et_noWalkers04
load et_noWalkers05
load et_noWalkers06
load et_noWalkers07
load et_noWalkers08
load et_noWalkers09 %et,nowalkers,e0

n_steps = 230000;

e_min = -102.915340231900; %h2@c60
% e_min = -27.5893786053000; %c12h10o
%e_min = -0.440620475522069; %peroskite
%e_min = -0.670556596200000; %h2
% e_min = -102.2448840269; %c60
% e_min = -5.69922832840000; %c2h6

start = 200000;

et_noWalkers_all = zeros(n_steps,3);
et_noWalkers_all(:,2) = et_noWalkers00(1:n_steps,2) + et_noWalkers01(1:n_steps,2) + et_noWalkers02(1:n_steps,2) + et_noWalkers03(1:n_steps,2) + et_noWalkers04(1:n_steps,2) + et_noWalkers05(1:n_steps,2) + et_noWalkers06(1:n_steps,2) + et_noWalkers07(1:n_steps,2) + et_noWalkers08(1:n_steps,2) + et_noWalkers09(1:n_steps,2);
et_noWalkers_all(:,3) = et_noWalkers00(1:n_steps,3) + et_noWalkers01(1:n_steps,3) + et_noWalkers02(1:n_steps,3) + et_noWalkers03(1:n_steps,3) + et_noWalkers04(1:n_steps,3) + et_noWalkers05(1:n_steps,3) + et_noWalkers06(1:n_steps,3) + et_noWalkers07(1:n_steps,3) + et_noWalkers08(1:n_steps,3) + et_noWalkers09(1:n_steps,3);
et_noWalkers_all(:,3) = et_noWalkers_all(:,3)/10;

all = zeros(10,n_steps,3);
all(1,:,:) = et_noWalkers00(1:n_steps,:);
all(2,:,:) = et_noWalkers01(1:n_steps,:);
all(3,:,:) = et_noWalkers02(1:n_steps,:);
all(4,:,:) = et_noWalkers03(1:n_steps,:);
all(5,:,:) = et_noWalkers04(1:n_steps,:);
all(6,:,:) = et_noWalkers05(1:n_steps,:);
all(7,:,:) = et_noWalkers06(1:n_steps,:);
all(8,:,:) = et_noWalkers07(1:n_steps,:);
all(9,:,:) = et_noWalkers08(1:n_steps,:);
all(10,:,:) = et_noWalkers09(1:n_steps,:);

means = zeros(10,1);
for i=1:10
    means(i) = mean(all(i,start:end,3))-e_min;
end 
ZPE = mean(means)
error = 0;
for i=1:10
    error = error + (means(i) - ZPE)^2;
end

error = sqrt(error/10)

st = std(et_noWalkers_all(start:end,3)-e_min)

start = 200;

figure(1)
x = 1:n_steps;
hold on
plot(x(start:end),(et_noWalkers00(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers01(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers02(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers03(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers04(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers05(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers06(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers07(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers08(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers09(start:n_steps,3)-e_min),'LineWidth',0.2)
plot(x(start:end),(et_noWalkers_all(start:n_steps,3)-e_min),'k','LineWidth',2.0)
xlabel('Iterations','Interpreter','latex')
ylabel('ZPE [Hartee]','Interpreter','latex')
grid on
xlim([0 230000])
% ylim([0.0725 0.076])
hold off

figure(2)
plot(x(start:end),et_noWalkers_all(start:end,2),'LineWidth',1.2)
xlabel('Iterations')
ylabel('Number of walkers')

figure(3)
plot(x(start:end),(et_noWalkers_all(start:end,3)-e_min),'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('$E_T - E_{min}$ [Hartee]','Interpreter','latex')


start = 2000;
le = n_steps - start;
err = zeros(le,1);
for i=1:le
    err(i) = abs(mean(et_noWalkers_all(i:i+start,3))-e_min-ZPE);
end

figure(4)
plot(1:le,err,'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('$Mean(E_T) - E_{min}$ [Hartee]','Interpreter','latex')
