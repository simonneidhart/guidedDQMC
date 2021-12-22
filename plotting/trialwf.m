load pes.out
load walker_positions.out
load new_masses.out
load eigenvalues.out
Bohr_Ang=0.5291772109;

dim = 180; %h2@c60
e_min = -102.915340231900; %h2@c60

e = reshape(pes,100,dim);
no_walkers = floor(length(walker_positions)/dim);
pos = reshape(walker_positions,dim,[])';

b = 25; %dont plot the first b and last b energies of the PES
sigma2 = zeros(dim,1);
x = ((1:100) - 51)*0.1/Bohr_Ang;
x_ = ((1:1000) - 510)*0.01/Bohr_Ang;

if (1)
    for i=1:5 %plot dimension 1 to 5
        figure(i)
        sigma2(i) = 1/sqrt(eigenvalues(i)*new_masses(i));
        y = 0.5*eigenvalues(i)*x.^2;
        y_gauss = gauss(x,sigma2(i));
        
        hold on
        histogram(pos(:,i),20,'Normalization','pdf')
        plot(x(b:end-b+1),e(b:end-b+1,i) - e_min,'b',x(b:end-b+1),y(b:end-b+1),'r',x(b:end-b+1),y_gauss(b:end-b+1),'g','LineWidth',1.1)
        legend('Walkers','PES','Harmonic approximation','Guiding wave function','Interpreter','latex')
        ylim([0 0.6])
        xlim([-3 3])
        grid('on')
        xlabel('$x$ [Angstrom]','Interpreter','latex')
        ylabel('Energy [Hartree]','Interpreter','latex')
        hold off

    end
end

function y = gauss(x,s2)
    y = 1/sqrt(2*pi*s2)*exp(-(x).^2/(2*s2));
end 