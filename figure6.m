%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is the code for figure 6 in Secure Quantum Remote Sensing without
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
precisionRoot = 10;
loopsRoot = 10;

%data extraction & re-organisation into appropriate array
fileID = fopen([userpath,sprintf('/nuMean_L%dP%d.bin',loopsRoot,precisionRoot)]);
nuMean = fread(fileID,'double');
fclose(fileID);
[sizeMatrix,~] = size(nuMean);
nuMean = reshape(nuMean,[4,10]);

%creation of grid for likelihood functions with 0 point & no 2pi point
precision = 2^precisionRoot;
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];
%plot variable choices. plotNumbers correspond to those done on HPC.
%Choices are indicative but with suitable time this can be run for many
%choices of plotNumber but change must also be made when creating data with
%HPC file.
plotNumbers = [30,60,100,150];
[~,plotNumbersNumber] = size(plotNumbers);
maxPasses = 10;

%figure creation
figure
hold on
for a = 1:plotNumbersNumber
    [~,xlinePoint] = min(nuMean(a,:));
    color_order = get(gca,'colororder');
    if a==1
        plot(1:maxPasses,nuMean(a,1:maxPasses),'-','color',color_order(a+2,:),'DisplayName',[num2str(plotNumbers(a)) , ' photons'],'Linewidth',3)
        xline(xlinePoint,'-','color',color_order(a+2,:),'HandleVisibility','off','Linewidth',1.5)
    elseif a==2
        plot(1:maxPasses,nuMean(a,1:maxPasses),'--','color',color_order(a+2,:),'DisplayName',[num2str(plotNumbers(a)) , ' photons'],'Linewidth',3)
        xline(xlinePoint,'--','color',color_order(a+2,:),'HandleVisibility','off','Linewidth',1.5)
    elseif a==3
        plot(1:maxPasses,nuMean(a,1:maxPasses),'-.','color',color_order(2,:),'DisplayName',[num2str(plotNumbers(a)) , ' photons'],'Linewidth',3)
        xline(xlinePoint,'-.','color',color_order(2,:),'HandleVisibility','off','Linewidth',1.5)
    elseif a==4
        plot(1:maxPasses,nuMean(a,1:maxPasses),'.-.','color',color_order(a+2,:),'DisplayName',[num2str(plotNumbers(a)) , ' photons'],'Linewidth',3)
        xline(xlinePoint,'.-.','color',color_order(a+2,:),'HandleVisibility','off','Linewidth',1.5)
    elseif a==5
        plot(1:maxPasses,nuMean(a,1:maxPasses),':','color',color_order(a+2,:),'DisplayName',[num2str(plotNumbers(a)) , ' photons'],'Linewidth',3)
        xline(xlinePoint,':','color',color_order(a+2,:),'HandleVisibility','off','Linewidth',1.5)
    else
    end   
end
hold off
ax = gca;
ax.FontSize = 18;
xlim([1 10])
xlabel('Number of passes','FontSize',18)
ylabel('Circular standard deviation','FontSize',18)
legend('Location','northwest','FontSize',18)

