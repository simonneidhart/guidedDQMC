load ener

ener = ener + 102.915322308400;

figure(1)
%x = (1:300)*0.01 - 0.01;
x = (1:100)*0.01 - 0.37189538;
plot(x,ener)
xlabel('rotation z axix [degree]')
%xlabel('shift in x direction [ang]')
ylabel('energy [hartee]')