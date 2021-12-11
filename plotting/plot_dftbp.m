load et_noWalkers_e0.out
load et_noWalkers
load walker_positions_kartesian

% et_noWalkers = et_noWalkers_e0;

e_min = -102.915340231900; %h2@c60
% e_min = -27.5893786053000; %c12h10o
%e_min = -0.440620475522069; %peroskite
%e_min = -0.670556596200000; %h2
% e_min = -102.2448840269; %c60
% e_min = -5.69922832840000; %c2h6
%unguided_res = 0.07386; 

figure(1)
ll = length(et_noWalkers);
x = 1:ll;
s = 200;
plot(x(s:end),(et_noWalkers(s:end,1)-e_min),'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('ZPE [Hartee]','Interpreter','latex')
grid('on')
% title('guiding wf adjusted @ n=800')

s = 200000;
zpe = mean(et_noWalkers(s:end,1))-e_min
st = std(et_noWalkers(s:end,1)-e_min)

figure(2)
plot(x,et_noWalkers(:,2),'LineWidth',1.2)
xlabel('Iterations')
ylabel('Number of walkers')

start = 2000;
le = ll - start;
err = zeros(le,1);

for i=1:le
    err(i) = abs(mean(et_noWalkers(i:i+start,3))-e_min-zpe);
end

figure(3)
plot(1:le,err,'LineWidth',1.2)
xlabel('Iterations','Interpreter','latex')
ylabel('$Mean(E_T) - E_{min}$ [Hartee]','Interpreter','latex')
grid('on')

nat = 62;
s = 10.0;

wp = walker_positions_kartesian;
walkers = et_noWalkers(end,2);

c = ones(nat,3);
c(end,:) = 2;
c(end -1,:) = 2;
c = repmat(c,walkers,1);

c_at = zeros(nat-2*walkers,3);
h_at = zeros(walkers,3);
h2_at = zeros(walkers,3);
c_count = 1;
h_count = 1;
h2_count = 1;
for i=1:nat*walkers
    if mod(i,62) == 0
        h_at(h_count,:) = wp(i-1,:);
        h_count = h_count + 1;
        h2_at(h2_count,:) = wp(i,:);
        h2_count = h2_count + 1;
    else
        c_at(c_count,:) = wp(i,:);
        c_count = c_count + 1;
    end 
end 

% h_at = wp(60:60:end,:);
% h2_at = wp(62:62:end,:);
        
       
figure(4)
hold on
scatter3(c_at(:,1),c_at(:,2),c_at(:,3),s,'b','filled')
% scatter3(wp(:,1),wp(:,2),wp(:,3),s,'b','filled')
scatter3(h_at(:,1),h_at(:,2),h_at(:,3),s,'y','filled')
scatter3(h2_at(:,1),h2_at(:,2),h2_at(:,3),s,'g','filled')
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
xlabel('x [Angstrom]','Interpreter','latex')
ylabel('y [Angstrom]','Interpreter','latex')
zlabel('z [Angstrom]','Interpreter','latex')
grid on
hold off

    
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





