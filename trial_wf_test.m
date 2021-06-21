load ener_hess
load eigenvalues
load walker_positions
load approx
load n_masses

%eigenvalues = eigenvalues_mass;
new_masses = n_masses;
e = reshape(ener_hess,100,18);
et = ones(1,100)*(-5.6253  + 5.6992283284);
no_walkers = floor(length(walker_positions)/18);
pos = reshape(walker_positions,18,no_walkers)';

b = 5;
kpot = zeros(18,1);
x = ((1:100) - 51)*0.1;

% for i=1:18
%     for j=1:100
%         if (e(j,i) < -5.6253)
%             kpot(i) = 2*(5.6992+e(j,i))/(abs(x(j)).^2);
%             break
%         end
%     end
% end

energy = 0;
for i=1:18
    %kpot(i) = 2*approx(i,2)/(approx(i,1))^2;
    %kpot(i) = sqrt(sqrt(eigenvalues(i)));
    kpot(i) = eigenvalues(i);
    energy = energy + 1/2*sqrt(kpot(i)/new_masses(i));
end 
energy

if (1)
    for i=1:18
    %     fit2 = fit(x(b:end-b+1)',e(b:end-b+1,i),'poly2');
    %     coeff = coeffvalues(fit2);
    %     kpot(i) = coeff(1);

        figure(i)
        y = 0.5*eigenvalues(i)*x.^2;
        y_new = 0.5*kpot(i)*x.^2;
        y_gauss = gauss(x,1/sqrt(kpot(i)*new_masses(i)));
    %     plot(x(b:end-b+1),e(b:end-b+1,i),x(b:end-b+1),y(b:end-b+1),x(b:end-b+1),y_new(b:end-b+1),x(b:end-b+1),et(b:end-b+1),'LineWidth',1.1)

        hold on
        histogram(pos(:,i),10,'Normalization','pdf')
        plot(x(b:end-b+1),e(b:end-b+1,i) + 5.6992283284,x(b:end-b+1),y_new(b:end-b+1),x(b:end-b+1),y_gauss(b:end-b+1),x(b:end-b+1),et(b:end-b+1),'LineWidth',1.1)
        legend('hist','actual','harmonic approx','trial wf','et')
        %legend('actual','harmonic approx','trial wf','et')
        ylim([0 1])
        xlabel('x')
        ylabel('Energy')
        hold off

    %     ylim([-5.7 -5.6])
    %     xlabel('x')
    %     ylabel('Energy')
    %     %plot(x(b:end-b+1),e(b:end-b+1,i),x(b:end-b+1),et(b:end-b+1))
    %     legend('actual','old approx','new approx','et')
    %     %legend('actual','et')
        title(int2str(i))
    end
end
save('k_pot','kpot','-ascii')

% figure(2)
% fitpoly2=fit(x(b:end-b+1)',e(b:end-b+1,18),'poly2')
% % Plot the fit with the plot method.
% plot(fitpoly2,x(b:end-b+1)',e(b:end-b+1,18))
% % Move the legend to the top left corner.
% legend('Location','NorthWest' );

function y = gauss(x,s2)
    y = 1/sqrt(2*pi*s2)*exp(-(x).^2/(2*s2));
end 