%% Load data
cd("~/Git/BART-RDD/Application/")
load("gpa.csv")
N = size(a.y);
N = N(1);
Y = a.y;
Z = a.w;
X = a.x;
c =0;
W = (X >= c);
%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50;
optTree.maxLevel        = 10;       % Maximum depth of the tree
optTree.maxLeaves       = 15;
optTree.minObs          = 50;       % Minimum observations within each leaf
optTree.cp              = 0.0001;   % Stopping rule for growing the tree: criterion must improve at least by this amount
optTree.maxIterGrow     = Inf;      % Maximum iteration during the growth of the tree
optTree.numKfold        = 5;        % number of K in K-fold cross-validation
optTree.CV1SE           = false;    % use of 1-SE rule in cross-validation
optTree.model           = 'linear'; % Set the model type: in case of RDD 'linear' is required
                                       % alternatively 'mean' can be used for rpart or Athey-Imbens type of
                                       % estimation, where there is only means computed
optTree.honest          = true;     % Set to use honest approach
optTree.type            = 'CATE';   % CATE is for sharp RD 'CLATE' id for fuzzy design
optTree.criteria        = 'RDD';    % Use of RDD criterion,
                                      % alternatively one can use 'athey-imbens' for experimental and observational studies with unconfoundedness
                                      % or 'MSE' for rpart prediction case
optTree.varEstimator    = 'hce-1';  % type of variance estimator: 'simple','hce-0','hce-1'.
                                      % For clustered SE, one needs to add the clusters when creating the sample
                                      % and use one of the estimators here as well
optTree.criteria_weight = 0.5;      % Criterion weight: 0.5 -> EMSE, using 1 -> `adaptive' criterion.
optTree.splitType       = 'causal-tree-RDD';
                                    % Set the splitting criterion to causal-tree RDD.
                                       % alternatively, one can use other type of splitting methods as well!
optTree.obsBucket       = 4;        % Minimum number of observations in the buckets
optTree.maxBucket       = Inf;      % Maximum number of buckets
optTree.orderPolinomial = 3;        % Order of polinomial used during the estimation
setCritFnc( optTree );              % Set everything in the object of optTree


%% Set-up the sample object
% Indexes for the training sample and for the estimation sample:
%   here just simply split the sample in the middle, but can use any
%   sampling tecnique.
if mod( N / 2 , 2 ) == 0
    indTrEst = [ 1 : N / 2 ; N / 2 + 1 : N ]';
else
    indTrEst = [ 1 : ( N + 1 ) / 2 ; [ ( N + 1 ) / 2 + 1 : N , N ] ]';
end
% Set up the sample
obj_sample = samples( optTree , Y , Z , 'index4TrEst' , indTrEst , 'X' , X , 'cutoff' , c );
% can use 'cluster' and the variable to estimate clustered SE

%% Find and estimate the optimal tree
[ est_tree , large_tree , opt_gamma , cv_crit , cand_gamma , cv_all ] = runTree( obj_sample , optTree );
% Outputs:
%  est_tree: estimated optimal RD tree
% not necessary outputs:
%  large_tree: large tree estimated on training sample
%  opt_gamma: optimal pruning parameter
%  cv_crit: averaged cross-validation criteria for each candidate penalty parameters
%  cand_gamma: candidate penalty parameters
%  cv_all: cross-validation criteria for each candidate penalty parameters and for all folds

%% Visualize the tree:
tostring( est_tree )
%% Predict tau
test = find(X>-0.1 & X<0.1);
Z_test = Z(test,:);
tau = predict(est_tree,Z_test);
resultfile = res_file;
writematrix(tau,resultfile)
