%%  Written by Sean William Moore 2022-11, updated for public use 2024-01-30, © CC 4.0

%This is computationally intensive code that produces the data for figures 3, 4 & 6 in Secure Quantum Remote Sensing without entanglement (2023) S.W.Moore, J.A.Dunnigham

%Please cite:
%   Sean W. Moore, Jacob A. Dunningham; Secure quantum remote sensing without entanglement. AVS Quantum Sci. 1 March 2023; 5 (1): 014406. https://doi.org/10.1116/5.0137260

%%  Introduction

%Variable choices available:
%   precisionRoot (multipass)
%   plotNumbers (multipass)
%   loopsRoot (multipass)
%   number of qubits (limited data)
%   loopsRoot (limited data)
%   plotPointsRoot (limited data)
%   precisionRoot (limited data)

%Will create files in same directory, otherwise use the following:
%userpath('/Users/you/folders')

%Reset random number generator 
%Add random number here before running on HPC, if running on local computer
%use sum(100*clock) 
reset(RandStream.getGlobalStream,sum(100*clock));

%%  Multipass

%basis grid creation for range [0,2*pi)
precisionRoot = 10;
precision = 2^precisionRoot;
%redeclared inside loops
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];

%shift on initial states from |X+>
chi = 0;
%shift in Bob's measurement basis from |X+>, pi/2 is Y basis.
eta = pi/2;

%Number of qubits
plotNumbers = [30,60,100,150];
[~,plotNumbersNo] = size(plotNumbers);

%true values over which we average likelihoods
plotPoints = 1:1:precision;
[~,plotPointsNo] = size(plotPoints);

%number of repetitions averaged over
loopsRoot = 10;
loops = 2^loopsRoot;

%maximum number of passes calcualted over
maxPasses = 10;
nuMean = zeros(plotNumbersNo,maxPasses);


for a=1:plotNumbersNo
    disp(['a= ',num2str(a)])
    plotPoints = 1:1:precision;
    for b=1:maxPasses %Number of passes in multipass test
        %to track progress if taking a long time
        disp(['a= ',num2str(a) , ' & b= ',num2str(b) ])

        %array declaration
        nu = zeros(1,plotPointsNo);

        %parallel processing 
        parfor c=1:plotPointsNo %Number of points in the plot to average over
            %array declaration
            likelihood = zeros(1,precision);
            %redeclaration to avoid overhead in parallel processing
            phi = linspace(0,2*pi,precision+1);
            phi(end) = [];

            for d=1:loops %Number of repetitions to average over
                %likelihood creation
                llh1 = FourStateLLHrnd(precision,phi(plotPoints(c)),chi,eta,plotNumbers(a),1);
                llhN = FourStateLLHrnd(precision,phi(plotPoints(c)),chi,eta,plotNumbers(a),b);
                %likelihood combination & renormalisation
                llh1N = llh1.*llhN;
                llh1N = llh1N/(sum(sort(llh1N*dphi)));
                %adding result to average
                likelihood = likelihood + llh1N/loops;
            end
            %finding circular standard deviation
            [~,nu(c)] = circStatistics(phi,likelihood);
        end
        %putting result in external array
        nuMean(a,b) = mean(nu);
    end
end
%exporting data to file
fileID = fopen([userpath,sprintf('/nuMean_L%dP%d.bin',loopsRoot,precisionRoot)],'w');
fwrite(fileID,nuMean,'double');
fclose(fileID);

%%  Circular standard deviation and Bias

%basis grid creation for range [0,2*pi)
precisionRoot = 15;
precision = 2^precisionRoot;
phi = linspace(0,2*pi,precision+1);
dphi = phi(2)-phi(1);
phi(end) = [];

%true values which will be used for plotting & array
plotPointsRoot = 8;
plotPointsDistance = precision/(2^plotPointsRoot);
plotPoints = 1:plotPointsDistance:precision;
[~,plotPointsNo] = size(plotPoints);
phi2 = linspace(0,2*pi,plotPointsNo+1);
phi2(end) = [];

%number of loops
loopsRoot = 15;
loops = 2^loopsRoot;

%shift on initial states from |X+>
chi = 0;
%shift in Bob's measurement basis from |X+>, pi/2 is Y basis.
eta = pi/2;

%number of photons calculated for
%n = [50,200,1000];
[~,nCount] = size(n);

%array declaration
nuMeanMatrix = zeros(nCount,plotPointsNo);

MLEMeanBiasMatrix = zeros(nCount,plotPointsNo);
for a=1:nCount
    %to track progress if taking a long time
    disp(['a= ',num2str(a)])

    %parallel processing
    parfor b=1:plotPointsNo
        %array declaration
        MLE = zeros(1,loops);
        nu = zeros(1,loops);
        
        %declaration within parfor to avoid overhead
        n = [50,200,1000];

        
        for c=1:loops
            %likelihood function creation
            likelihood = FourStateLLHrnd(precision,phi2(b),chi,eta,n(a));
            %find circular standard deviation
            [~,nu(c)] = circStatistics(phi,likelihood);
            %find MLE & record
            [~,maxIndex] = max(likelihood);
            MLE(c) = phi(maxIndex);
        end

        %putting result in external array
        nuMeanMatrix(a,b) = mean(nu);

        %finding MLE mean & attributing difference with true value to external array
        MLEMean = circStatistics(MLE);
        MLEMeanBiasMatrix(a,b) = wrapToPi(MLEMean - phi(plotPoints(b)));
    end
end

%Saving data to files
fileID = fopen([userpath,sprintf('/nuMeanMatrix_L%dN%dPP%d.bin',loopsRoot,plotPointsRoot,precisionRoot)],'w');
fwrite(fileID,nuMeanMatrix,'double');
fclose(fileID);

fileID = fopen([userpath,sprintf('/MLEMeanBiasMatrix_L%dN%dPP%d.bin',loopsRoot,plotPointsRoot,precisionRoot)],'w');
fwrite(fileID,MLEMeanBiasMatrix,'double');
fclose(fileID);




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








%%  Circular statistics

function [mu , nu, var ] = circStatistics(alpha, w)
%individual data points
%alpha      ~data point
%w          ~weight of data point
%histograms:
%alpha      ~basis of histogram
%w          ~values on histogram

%output:
%mu     ~mean direction
%nu     ~circular standard deviation
%var    ~circular variance

if nargin < 2 || isempty(w) 
%assume single data points of equal weight
	w = ones(size(alpha));
end

%Find mean resultant vector
z = sum(w.*exp(1i*alpha));

%obtain mean by direction in range [0,2*pi)
mu = mod(angle(z),2*pi);

%mean resultant length
Rbar = abs(z)/sum(w);

%circular standard deviation
nu = sqrt(-2*log(Rbar));

%circular variance, two competing definitions:

%More common in directional statistics: 0 no variance, 1 maximal variance
%var = 1 - Rbar;

%Tends to variance for large data/low dispersion
var = 2*(1 - Rbar);

end

%

%Avoid instances where numerical noise causes sum(p)~=1
function r = mnrnd2(n,p)
if sum(p)~=1
    p = p./sum(p);
end
r = mnrnd(n,p);
end
