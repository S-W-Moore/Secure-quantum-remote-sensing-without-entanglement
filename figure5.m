%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is the code for figure 5 in Secure Quantum Remote Sensing without
%entanglement (2023) S.W.Moore, J.A.Dunnigham

%Please cite:
%   Sean W. Moore, Jacob A. Dunningham; Secure quantum remote sensing without entanglement. AVS Quantum Sci. 1 March 2023; 5 (1): 014406. https://doi.org/10.1116/5.0137260

%reset random numbers
reset(RandStream.getGlobalStream,sum(100*clock)); 
color_order = get(gca,'colororder');  
%%  Multipass method 30,4

%basis grid creation for range [0,2*pi)
precisionRoot = 10;
precision = 2^precisionRoot;
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];

%choice of true parameter value, can be made random in range [0,2*pi)
phi0 = 1.25;%phi(phiR);

%shift on initial states from |X+>
chi = 0;
%shift in Bob's measurement basis from |X+>
eta = pi/2;

%number of qubits in single pass & multipass to be combined. Half the total
Ntot = 30;

%array declaration
llh1 = zeros(1,precision);
llhN = zeros(1,precision);
llh1N = zeros(1,precision);

llh1doub = zeros(1,precision); %For putting all photons through single pass

%number of passes
N=4;

%summing of many possible likelihoods
loops = 1000;
parfor a=1:loops
%1 pass
    l1 = FourStateLLHrnd(precision,phi0,chi,eta,Ntot,1);
%N passes
    lN = FourStateLLHrnd(precision,phi0,chi,eta,Ntot,N);
%combination    
    l1N = l1.*lN;
%renormalisation
    l1N = l1N/(sum(sort(l1N)*dphi));
%addition to sum of likelihoods
    llh1 = llh1 + l1;
    llhN = llhN + lN;
    llh1N = llh1N + l1N;
%double data single pass, add to sum of likelihoods
    llh1doub = llh1doub + FourStateLLHrnd(precision,phi0,chi,eta,2*Ntot,1);
end

%take the averages from the sums
llh1 = llh1/loops;
llhN = llhN/loops;
llh1N = llh1N/loops;
llh1doub = llh1doub/loops;
    
%plotting
figure
plot(phi/pi,llh1,'--','color',color_order(3,:),'Linewidth',3,'DisplayName','1 pass 30 photons')
hold on
plot(phi/pi,llhN,':','color',color_order(6,:),'Linewidth',3,'DisplayName',[num2str(N) , ' passes 30 photons'])
plot(phi/pi,llh1N,'color',color_order(4,:),'Linewidth',3,'DisplayName',['Combined 1&',num2str(N), ' passes'])
plot(phi/pi,llh1doub,'-.','color',color_order(2,:),'Linewidth',3,'DisplayName','1 pass 60 photons')
xline(phi0/pi,'Linewidth',1.5,'DisplayName','True Value')
hold off
ax = gca;
ax.FontSize = 18;
ylabel('Likelihood function','FontSize',18)
xlabel('\phi/\pi','FontSize',18)
legend('FontSize',18)



%   This function uses a multinomial to make produce data for SQRSWE and
%   outputs a grid approximation to the likelihood function. It supports
%   more than 1000 qubits & multiple passes.
function [likelihood, N] = FourStateLLHrnd(precision,phi0,chi,eta,Ntot,passes)
%   precision = number of data points for the produced likelihood function
%   phi0 = true value
%   chi = angle that the states arriving at Bob are offset by compared to |X+>
%   eta = angle of Bob's measurement compared the |X+>
%   Ntot = number of photons entering the system
%   passes = number of passes
if nargin<6
    passes = 1;
end

%probability array creation
P = [ (1+ cos(passes*phi0+chi-eta)) , (1+ cos(passes*phi0+chi-eta+pi/2)) , (1+ cos(passes*phi0+chi-eta+pi)) , (1+ cos(passes*phi0+chi-eta+3*pi/2)) ]/4;

%Result array creation
N =  mnrnd2(Ntot,P);
if any(isnan(N), 'all')
    error('NaN errors!')
end

%Adaptation to high N. For large N double precision arithmatic is
%insuffucuent causing NaN/Inf errors without this step
Nmultiple = 1;
if Ntot>1000
    Nmultiple = Ntot/1000;
    N = N /Nmultiple;
end
    
%recreation of grid
dphi = 2*pi/precision;
phi = linspace(0,2*pi,precision+1);
phi(end) = [];

%Likelihood function grid approximaiton creation
L1 = .5*(1+ cos(passes*phi+chi-eta));
L2 = .5*(1+ cos(passes*phi+chi-eta+pi/2));
L3 = 1-L1;
L4 = 1-L2;
likelihood = L1.^N(1).*L2.^N(2).*L3.^N(3).*L4.^N(4);
likelihood = likelihood/(sum(sort(likelihood))*dphi);
if Ntot>1000
    likelihood = likelihood.^Nmultiple;
    likelihood = likelihood/(sum(sort(likelihood))*dphi);
end

end

%Avoid instances where numerical noise causes sum(p)~=1
function r = mnrnd2(n,p)
if sum(p)~=1
    p = p./sum(p);
end
r = mnrnd(n,p);
end