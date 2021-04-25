load ener
load geom_opt.x
load walker_positions_unguided
e = reshape(ener,100,24);
et = ones(1,100)*(-5.6253  + 5.6992283284);
no_walkers = length(walker_positions_unguided)/8;
pos = reshape(walker_positions_unguided',24,no_walkers)';

masses = [12.011, 12.011, 12.011, 12.011, 12.011, 12.011, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007];
masses = masses.*1836.1;

b = 5;
kpot = zeros(24,1);
x = ((1:100) - 51)*0.1;

% for i=1:18
%     for j=1:100
%         if (e(j,i) < -5.6253)
%             kpot(i) = 2*(5.6992+e(j,i))/(abs(x(j)).^2);
%             break
%         end
%     end
% end

% for i=1:18
%     kpot(i) = 2*approx(i,2)/(approx(i,1))^2;
% end 

st = zeros(24,1);
for i=1:24
%     fit2 = fit(x(b:end-b+1)',e(b:end-b+1,i),'poly2');
%     coeff = coeffvalues(fit2);
%     kpot(i) = coeff(1);
    
%     y = 0.5*eigenvalues(i)*x.^2;
%     y_new = 0.5*kpot(i)*x.^2;
%     y_gauss = gauss(x,1/sqrt(kpot(i)));
%     plot(x(b:end-b+1),e(b:end-b+1,i),x(b:end-b+1),y(b:end-b+1),x(b:end-b+1),y_new(b:end-b+1),x(b:end-b+1),et(b:end-b+1),'LineWidth',1.1)
    
%     hold on
%     histogram((pos(:,i)-ones(length(pos(:,i)),1)*inital(i)),10,'Normalization','pdf')
%     plot(x(b:end-b+1),e(b:end-b+1,i) + 5.6992283284,x(b:end-b+1),y_new(b:end-b+1),x(b:end-b+1),y_gauss(b:end-b+1),x(b:end-b+1),et(b:end-b+1),'LineWidth',1.1)
%     legend('hist','actual','harmonic approx','trial wf','et')
%     ylim([0 1])
%     xlabel('x')
%     ylabel('Energy')
%     hold off
    me = mean(pos(:,i));
    st(i) = var(pos(:,i));
    y_trial = gauss(x,st(i));
    y_approx = 0.5/st(i)^2*x.^2;
    
    
    figure(i)
    hold on
    xlabel('x')
    ylabel('Energy')
    histogram((pos(:,i)-ones(length(pos(:,i)),1)*me),20,'Normalization','pdf')
    plot(x(b:end-b+1),e(b:end-b+1,i)+ 5.6992283284 ,x(b:end-b+1),et(b:end-b+1),x(b:end-b+1),y_trial(b:end-b+1),x(b:end-b+1),y_approx(b:end-b+1))
    %legend('actual','old approx','new approx','et')
    legend('hist','epot','et','trial wf','harmonic approx')
    ylim([0 1])
    title(int2str(i))
end

st = reshape(st,3,8)';

save('sigma','st','-ascii')

% figure(2)
% fitpoly2=fit(x(b:end-b+1)',e(b:end-b+1,18),'poly2')
% % Plot the fit with the plot method.
% plot(fitpoly2,x(b:end-b+1)',e(b:end-b+1,18))
% % Move the legend to the top left corner.
% legend('Location','NorthWest' );

function y = gauss(x,s2)
    y = 1/sqrt(2*pi*s2)*exp(-(x).^2/(2*s2));
end 