%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is the code for figure 7 in Secure Quantum Remote Sensing without
%entanglement (2023) S.W.Moore, J.A.Dunnigham

%Please cite:
%   Sean W. Moore, Jacob A. Dunningham; Secure quantum remote sensing without entanglement. AVS Quantum Sci. 1 March 2023; 5 (1): 014406. https://doi.org/10.1116/5.0137260

reset(RandStream.getGlobalStream,sum(100*clock)); 
color_order = get(gca,'colororder');   

%%  Comparison of different cryptographic methodogies security against photon splitting attacks

%grid creation
x = 0:0.001:1;
%BB84
BB84 = (exp(x)-1 -x)./(exp(x)-1);
%secure quantum remote sensing
YQ4 = (exp(x)-1 -x -3*x.^2./8)./(exp(x)-1);

%plotting
figure
plot(x,BB84,'DisplayName','BB84','Linewidth',3,'color',color_order(4,:))
hold on
plot(x,YQ4,'--','DisplayName','SQRS','Linewidth',3,'color',color_order(6,:))
hold off
ax = gca;
ax.FontSize = 18;
xlabel('Mean photon number','FontSize',18)
ylabel('Relative information rate','FontSize',18)
legend('location','northwest','FontSize',18)