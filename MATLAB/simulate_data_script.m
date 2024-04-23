% Summary statistics for testing burdenEM, based on simulations in Weiner &
% Nadig et al

clear

% Change if needed
addpath('../../bhr/MATLAB/')

% Which replicate number (used for naming sumstats files)
rep = 1;

rng(rep, 'twister') % for reproducibility

% This script will produce one gene sumstats file for each simulation,
% where the number of simulations is 
% length(h2BurdenTrue) * length(max_af) * length(titles)
save_path = '../data/';
if ~isfolder(save_path)
    mkdir(save_path)
end

gg = 1.8e4; % no. genes [identical to Weiner & Nadig et al]

% number of traits under pleiotropic selection [identical to Weiner & Nadig
% et al]
noTraitsUnderSelection = 100;

% no. snps per gene [identical to Weiner & Nadig et al]
mm_per_gene = randi(1e3,gg,1);

% true h2 burden; more non-null values
h2BurdenTrue = [0 0.005 0.01];

% AF bins; unlike Weiner & Nadig et al, just one bin
max_af = 1e-4; % maximum MAF
min_af = 0; % minimum AF

% Width of effect-size distribution (2 narrow, 10 wide)
width = 5 * ones(1,10);

% Description of each simulation. Parameters below correspond to each
% simulation. [Fewer than W&N]
titles = {'realistic','small_N','large_N','strong_selection',...
    'strong_popstrat','more_polygenic','less_polygenic'};

% sample size
nn = [5e5 1e5 2e6 5e5 5e5 5e5 5e5]; 

% mean strength of selection (Ns)
selectionStrength = 10 * [1 1 1 10 1 1 1];

% Median of gene heritability distribution
sparsity = [2e3 * ones(1,5) 5e2 1e4];

% overdispersion effects
overdisp = zeros(size(nn));

% population stratification effects
popstratmean = [1 1 1 1 10 1 1] * 1e-5;
popstratvar = [1 1 1 1 10 1 1] * 1e-7;


clear h2est_noAnnot h2est_Annot h2est_AnnotObs
for nval=1:length(nn) % loop over simulation setups (titles)
    
    % mixture component variances
    sigmasq = [1/width(nval) 1 width(nval) 0];
    
    % probability each mixture cpt (cols) in each annotation (rows), with sum
    % at most 1
    prior = 1/sparsity(nval) ./ sigmasq(1:end-1);
    if sum(prior)>1; error('prior should sum to <=1'); end
    prior = [prior, 1 - sum(prior)]; % null component
    
    for sim = 1:length(h2BurdenTrue) % loop over true h2 values
        for bin = 1:length(max_af) % loop over AF bins
            
            % simulate gene/variant level sumstats
            [genes, variants] = ...
                simulate_rare_sumstats(nn(nval),gg,mm_per_gene,sigmasq,...
                'overdispSupport',overdisp(nval)*sigmasq,...
                'sigmasqPrior',prior,...
                'maxAF',max_af(bin),...
                'minAF',min_af(bin),...
                'meanNs', 0,...
                'popStratVar',popstratvar(nval),...
                'popStratMean',popstratmean(nval),...
                'h2Target',h2BurdenTrue(sim),...
                'meanNs', selectionStrength(nval),...
                'noTraits',noTraitsUnderSelection,...
                'selectionModel','stabilizing');
            
            % save to file
            writetable(struct2table(genes),[save_path,titles{nval},'.MAFbin=',num2str(bin),'.h2=',num2str(h2BurdenTrue(sim)),'.rep=',num2str(rep),'.csv']);
            
        end
    end
    
end

