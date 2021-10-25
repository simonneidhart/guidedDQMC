load walker_positions_hessian_trial
Bohr_Ang=0.529177;
dim = 18;
now = 869;

pos = reshape(walker_positions_hessian_trial,dim,[])';

m = zeros(dim,1);
s2 = zeros(dim,1);
for i=1:dim
    m(i) = mean(pos(:,i));
    s2(i) = var(pos(:,i));
end 

x = ((1:100) - 51)*0.01/Bohr_Ang;
for i=1:18
    figure(i)
    y_gauss = gauss(x,s2(i));
    hold on
    histogram(pos(:,i),10,'Normalization','pdf')
    plot(x,y_gauss,'LineWidth',1.1)
    ylim([0 5])
    xlabel('x [Angstroem]','Interpreter','latex')
    ylabel('Energy [Hartree]','Interpreter','latex')
    hold off
end




function y = gauss(x,s2)
    y = 1/sqrt(2*pi*s2)*exp(-(x).^2/(2*s2));
end 