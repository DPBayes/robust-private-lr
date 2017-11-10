function [params] = initial_params()
    params = struct();
    params.synth_data=0; % 1 to use simulation data or 0 to use sanger data
    
    
    params.codir='../corr/'; % create the folder 
    params.wodir='../wpc/';
    params.todir='../tensor/';
    
    params.processed_data_filename = 'all_data.mat';%preprocessed data filename
    
    params.features=8;%bit setting flag for different genomic features GE = 4; CN = 3; MT = 2; TT = 1;
    params.corr_drug = 0;%corr_drug:0 ranking patients; 1 ranking drugs
    params.debug=0;
    params.sdebug=0;
    %feature selection method
    params.dimension_reduce =4;%dimension reduction method
    params.GE_feature_file='lasso_features_114.mat';
        
    params.real2rank=0;
    params.rmse=1;
    params.corr=1;
    params.cindex=0;
    params.wpc=1; %1: wpc, 0:cindex
    params.r2=1;
    params.rank_pat=1; %1 rank patients, 0 rank drugs 
    params.topK=5;
    params.chk_model=0;
    params.fold=50; %set as 5 for quick check, value used is 50 in paper
    
    if(params.synth_data)
        params.data_filename = 'synth_data.mat';%preprocessed data filename
        params.dimension_reduce=0;
    end
    params.no_models=6;
    
    % parameters
    params.sigma=1;
    params.theta=1;
    
    
  %  params.lambda=1e-2;%regularization value RR worse 
    params.lambda=1e0; 
  %  params.hs_sigma=1e-2;
  
    %params.sigma=1e-1;
    %params.theta=1e-1;
    params.lambda=1e0; 
    
    
    params.cross_valid=1;%need cross validataion or not
    if(params.debug || params.sdebug) params.cross_valid=0; end
    
    %dealing with missing values
    params.mod_nan=0; % 0: do not do anything for NaN explicitly
    params.pad=0; % replace NaN with 0
    
    %normalize
    params.logscaling = 0;%log transform and then scaling the feature vectors [0 better]
    
    params.l1normalize = 0;%l1-normalize [0 better]
    params.l1norm_dim = 2; % n x d matrix     
    
    params.l2normalize = 1;%l2-normalize [1 better]  
    params.l2norm_dim = 2; % n x d matrix     
    params.normalize=1;
    
    %test and train data ratio
    params.frac = 0.8;%the ratio for train data
    
    if(params.synth_data)
        params.total_size=1e3;
    else
        params.total_size=650; %ver1
        params.total_size=613; %ver2
        params.total_size=985; %ver3
    end
    params.tesize=round(params.total_size * (1 - params.frac)); %ver1
    params.tesize=100; %ver2,ver3
    
    %private data
    params.insize=10; %ver2,ver3 
    params.min_insize=0; %to set limit for T
    params.max_insize=30; %to fix final size of private data
    params.min_d=5;
    params.step_d=5;
    params.max_d=64; %ver3
    
    
    params.prsize = params.total_size - params.tesize - params.insize; %ver1
    
    %params.insize=13; %ver2 
    params.insize=10; %ver2, ver3
    params.prsize = params.total_size - params.tesize - params.max_insize; 
    %round(params.total_size*0.7); %-100-params.tesize; %500;%the number of private data points
    
    params.size_var = 1;%check the effect of the size of (private) data

    
   
    %privacy
    params.privacy_var = 1;%the default step is 0.2 for testing the trade-off between utility and privacy
    params.privacy_iter = 1;%generate multiple perturbation results to reduce the uncertainty of randomness noisy, by default:1
    params.noise_only=0;  
    
   % params.seed=10000;
    
    params.scale=1;
    params.mult_size=10; % working till 10
    params.extrapolate=0;
    if(params.extrapolate==0)
         params.mult_size=0;
    end
    params.intercept=0; %[0 better]
   
    params.no_internal=0;
    params.only_internal =0;
    params.denoise=0;
    params.ITER=1;
    
     params.old_normalizey=0;
    params.anal_data=0;
    
    params.only_plot=0;
    %params.odir=[params.odir,'OP/'];
    %params.odir=[params.odir,'P3/'];
    %params.odir=[params.odir,'FM/'];
    %params.odir=[params.odir,'OPR/'];
    %params.odir=[params.odir,'P/'];
    %params.odir=[params.odir,'L/'];
     % params.odir=[params.odir,'LP/'];
    %params.tensor_path=[params.odir,'corr/'];
end
    
  
