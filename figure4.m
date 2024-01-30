%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is the code for figure 4 in Secure Quantum Remote Sensing without
%entanglement (2023) S.W.Moore, J.A.Dunnigham

%Please cite:
%   Sean W. Moore, Jacob A. Dunningham; Secure quantum remote sensing without entanglement. AVS Quantum Sci. 1 March 2023; 5 (1): 014406. https://doi.org/10.1116/5.0137260

reset(RandStream.getGlobalStream,sum(100*clock)); 
color_order = get(gca,'colororder');  

%%  HPC data extraction
%Example data produced using a high performance cluster is avaiable in the
%folder, as is the code used to produce the data.

%update userpath to folder where file is download unless you are in the
%same workspace

%userpath('/Users/you/folders')

%set the following two variable to align with HPC data, different
%calculation precisions for different running times in HPC code. No need to
%change unless you have produced your own data.
precisionRoot = 15;
plotPointsRoot =8;
plotPointsNo = 2^plotPointsRoot;
loopsRoot = 15;

%data extraction & re-organisation into appropriate array
fileID = fopen([userpath,sprintf('/MLEMeanBiasMatrix_L%dN%dPP%d.bin',loopsRoot,plotPointsRoot,precisionRoot)]);
MLEMeanBiasMatrix = fread(fileID,'double');
fclose(fileID);
MLEMeanBiasMatrix = reshape(MLEMeanBiasMatrix,[3,plotPointsNo]);

%creation of grid for likelihood functions with 0 point & no 2pi point
phi = linspace(0,2*pi,plotPointsNo+1);
phi(end) = [];

%numbers corresponding to the amount of qubits used to produce the data
n=[50,200,1000];
%figure creation
figure
plot(phi(1:2:end)/pi,MLEMeanBiasMatrix(1,1:2:end),'--','color',color_order(3,:),'DisplayName',[num2str(n(1)), ' photons'],'Linewidth',3)
hold on
plot(phi(1:2:end)/pi,MLEMeanBiasMatrix(2,1:2:end),':','color',color_order(4,:),'DisplayName',[num2str(n(2)), ' photons'],'Linewidth',3)
plot(phi(1:2:end)/pi,MLEMeanBiasMatrix(3,1:2:end),'-','color',color_order(6,:),'DisplayName',[num2str(n(3)), ' photons'],'Linewidth',3)
hold off
xlabel('\phi/\pi','FontSize',18)
ylabel('Mean bias','FontSize',18)
ax = gca;
ax.FontSize = 18;
legend('location','northwest','FontSize',18)
ylim([-.07,0.07])