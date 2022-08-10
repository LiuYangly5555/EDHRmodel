function g = Graphdata(sim, day)

x=linspace(1,day+1,day*144);

subplot(3,1,1);
plot(x, sim.q1+sim.p1, x, sim.q2-sim.p1);
axis([1,day+1,-10, 20]);xticks(1:1:21);
grid on;xlabel("day");ylabel("Q");
legend('Q_{upper}', 'Q_{lower}', 'Location','northeast', 'Orientation', 'horizontal')

subplot(3,1,2);
plot(x, sim.SWC1, x, sim.SWC2);
axis([1,day+1,0,0.5]);xticks(1:1:21);
grid on;xlabel("day");ylabel("SWC");
legend('SWC_{shallow}', 'SWC_{deep}', 'Location','northeast', 'Orientation', 'horizontal')

subplot(3,1,3);
semilogy(x, sim.H_root_1,  x, sim.H_root_2, x, sim.H1,  x, sim.H2);
axis([1,day+1,-100000,-1]);xticks(1:1:21);
grid on;xlabel("day");ylabel("H");
legend('H_{root,upper}', 'H_{root,lower}', 'H_{soil,shallow}', 'H_{soil,deep}', 'Location','northeast', 'Orientation', 'horizontal')


end

