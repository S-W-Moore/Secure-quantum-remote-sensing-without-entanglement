%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is the code for figure 9 in Secure Quantum Remote Sensing without
%entanglement (2023) S.W.Moore, J.A.Dunnigham

%Please cite:
%   Sean W. Moore, Jacob A. Dunningham; Secure quantum remote sensing without entanglement. AVS Quantum Sci. 1 March 2023; 5 (1): 014406. https://doi.org/10.1116/5.0137260


reset(RandStream.getGlobalStream,sum(100*clock)); 
color_order = get(gca,'colororder');   
%%  MIM simulation

%Preparation of grid for likelihood function grid approximation
precisionRoot = 10;
precision = 2^precisionRoot;
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];

phi0 = .4*pi;   %change to 2*pi*rand for random true parameter
chi=0;          %Overall shift in state basis
eta=pi/2;       %Bob parameter estimation emasurement basis
passes = 1;     %Change for multipass 
NE = 1;         %Number of states that Eve manages to attack

%result array declaration
Allh20 = zeros(1,precision);
Allh90 = zeros(1,precision);
Allh380 = zeros(1,precision);
Ellh = zeros(1,precision);

%averaging over many calculations
loops = 100000;
parfor a=1:loops
    [temp1,~] = MIMattackLLHs(precision,phi0,chi,eta,passes,20,NE);
    Allh20 = Allh20 + temp1/loops;
    [temp2,~] = MIMattackLLHs(precision,phi0,chi,eta,passes,90,NE);
    Allh90 = Allh90 + temp2/loops;
    [temp3,temp4] = MIMattackLLHs(precision,phi0,chi,eta,passes,380,NE);
    Allh380 = Allh380 + temp3/loops;
    Ellh = Ellh + temp4/loops;
end

%plotting
figure
plot(phi/pi,Allh380,'-.','color',color_order(2,:),'DisplayName','Alice, p=0.95','Linewidth',3)
hold on
plot(phi/pi,Allh90,'--','color',color_order(3,:),'DisplayName','Alice, p=0.9','Linewidth',3)
plot(phi/pi,Allh20,'-','color',color_order(4,:),'DisplayName','Alice, p=0.8','Linewidth',3)
plot(phi/pi,Ellh,':','color',color_order(6,:),'DisplayName','Eve','Linewidth',3)
hold off
xline(0.4,'.-.','DisplayName','True value','Linewidth',1.5)
xlabel('\phi/\pi','FontSize',18)
ylabel('Mean likelihood','FontSize',18)
ax = gca;
ax.FontSize = 18;
legend('location','northeast','FontSize',18)


%When Eve attacks and hides the attack on some qubits during Alice's
%estimation process
function [Alice,Eve] = MIMattackLLHs(precision,phi0,chi,eta,passes,NA,NE)
%precision  ~grid size
%phi0       ~true value
%chi        ~overall shift in basis
%eta        ~basis of Bob's measurement
%eta        ~number of passes by Bob
%NA         ~total number of results for Alice
%NE         ~number of attacks Eve is hiding

if NE>NA
    error('MIMattackLLHs: NE>NA')
end

%How many times deos Eve choose to measure in the same basis as the initial
%state
NEnet = mnrnd(NE,[0.5,0.5]);
NEnet = NEnet(1);
%recreation of standard grid with 0 point & no 2pi point
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];

%Probability array for results
P = [ (1+cos(passes*phi0+chi-eta)) , (1+cos(passes*phi0+chi-eta+pi/2)) , (1+cos(passes*phi0+chi-eta+pi)) , (1+cos(passes*phi0+chi-eta+3*pi/2)) ]/4;

%Creation of Alice's trues results
N=  mnrnd2(NA-NEnet,P);
%Attribution of noise caused by Eve
if NE>0
    for b=1:NE
        temp = randi(4);
        N(temp) = N(temp) + 1;
    end
end

if any(isnan(N), 'all')
    error('NaN errors!')
end

%Creation of arrays for each probability possibility
L1 = .5*(1+cos(passes*phi+chi-eta));
L2 = .5*(1+cos(passes*phi+chi-eta+pi/2));
L3 = 1-L1;
L4 = 1-L2;

%Alice likelihood creation
likelihood = L1.^N(1).*L2.^N(2).*L3.^N(3).*L4.^N(4);
likelihood = likelihood/(sum(sort(likelihood))*dphi);

Alice = likelihood;

%Eve result creation
NEve =  mnrnd2(NE,P);
%Eve likelihood creation
likelihoodE = L1.^NEve(1).*L2.^NEve(2).*L3.^NEve(3).*L4.^NEve(4);
likelihoodE = likelihoodE/(sum(sort(likelihoodE))*dphi);

Eve=likelihoodE;

end

%Avoid instances where numerical noise causes sum(p)~=1
function r = mnrnd2(n,p)
if sum(p)~=1
    p = p./sum(p);
end
r = mnrnd(n,p);
end
