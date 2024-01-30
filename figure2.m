%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is the code for figure 3 in Secure Quantum Remote Sensing without
%entanglement (2023) S.W.Moore, J.A.Dunnigham

%Please cite:
%   Sean W. Moore, Jacob A. Dunningham; Secure quantum remote sensing without entanglement. AVS Quantum Sci. 1 March 2023; 5 (1): 014406. https://doi.org/10.1116/5.0137260

%reset random numbers
reset(RandStream.getGlobalStream,sum(100*clock)); 
color_order = get(gca,'colororder');  

%%  Plot of single true value mean likelihoods

%shift on initial states from |X+>
chi = 0;
%shift in Bob's measurement basis from |X+>, pi/2 is Y basis.
eta = pi/2;

%basis grid creation for range [0,2*pi)
precisionRoot = 12;
precision = 2^precisionRoot;
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];

%array declaration
llhA = zeros(1,precision);
llhS = zeros(1,precision);
llhC = zeros(1,precision);

%true parameter value
phi0 = 2.5;

%number of qubits
Ntot = 100;

%summing of many possible likelihoods
loops = 1000;
parfor a=1:loops
    llh = FourStateLLHrndType(precision,phi0,chi,eta,Ntot,1,'A');
    llhA = llhA + llh;
    llh = FourStateLLHrndType(precision,phi0,chi,eta,Ntot,1,'S');
    llhS = llhS + llh;
    llh = FourStateLLHrndType(precision,phi0,chi,eta,Ntot,1,'C');
    llhC = llhC + llh;
end

%taking average from sum
llhA = llhA/loops;
llhS = llhS/loops;
llhC = llhC/loops;

%plotting
figure
plot(phi/pi,llhA,'color',color_order(3,:),'DisplayName','\sigma_x & \sigma_y states','LineWidth',3)
hold on
plot(phi/pi,llhC,'-.','color',color_order(4,:),'DisplayName','\sigma_y states','LineWidth',3)
plot(phi/pi,llhS,'--','color',color_order(6,:),'DisplayName','\sigma_x states','LineWidth',3)
xline(phi0/pi,':','DisplayName','True value','LineWidth',3)
hold off
ax = gca;
ax.FontSize = 18;
xlabel('Parameter (\pi)')
ylabel('Likelihood function')
legend



%   This function uses a multinomial to make produce data for SQRSWE and
%   outputs a grid approximation to the likelihood function. It supports
%   more than 1000 qubits & multiple passes.

%   Type allows choice of using only cosine/sine or all types of
%   probabilities
function [likelihood, N] = FourStateLLHrndType(precision,phi0,chi,eta,Ntot,passes,type)
%   precision = number of data points for the produced likelihood function
%   phi0 = true value
%   chi = angle that the states arriving at Bob are offset by compared to |X+>
%   eta = angle of Bob's measurement compared the |X+>
%   Ntot = number of photons entering the system
%   passes = number of passes
if nargin<6
    passes = 1;
end
if nargin<6
    type = 'A';
end

%recreation of grid
dphi = 2*pi/precision;
phi = linspace(0,2*pi,precision+1);
phi(end) = [];

%probability array creation
if type == 'A' 
    P = [ .5*(1+ cos(passes*phi0+chi-eta)) , .5*(1+ cos( passes*phi0+chi-eta+pi/2)) , .5*(1+ cos( passes*phi0+chi-eta+pi)) , .5*(1+ cos( passes*phi0+chi-eta+3*pi/2)) ]/2;
    N =  mnrnd2(Ntot,P);
elseif type == 'S'
    P = [ .5*(1+ cos( passes*phi0+chi-eta)) , 0 , .5*(1+ cos( passes*phi0+chi-eta+pi)) , 0 ];
    N =  mnrnd2(Ntot,P);
elseif type == 'C'
    P = [ 0 , .5*(1+ cos( passes*phi0+chi-eta+pi/2)) , 0 , .5*(1+ cos( passes*phi0+chi-eta+3*pi/2)) ];
    N =  mnrnd2(Ntot,P);
else
    error( 'type value invalid')
end


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