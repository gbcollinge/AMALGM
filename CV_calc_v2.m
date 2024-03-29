function [ CV_score, std_dev] = CV_calc_v2( test_mat,energy,a,optional_reference_CV_score)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% CV_calc.m: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given an N x M submatrix corresponding to the M interaction coefficients
% of some (size M) cluster expansion (CE) for some N known DFT energies, 
% calculates the Leave Multiple Out - Cross Validation (LMO-CV) score.
%
% The training set of N DFT energies is split into a construction set (size
% N-d) and validation set (size d). The relative size of these is
% determined by the user specified "a": the fractional portion of the N
% total data that becomes the d validation set data.
% Per K. Baumann, Trends Anal. Chem. 22 (2003) 395-406, "a" should not be
% set below 0.4
%
% v2: Enforces an unbiased sampling of the data set. 
%     ...After N random cuts of the data, the algorithm chooses the
%     validation set from structures that have been biased against in the
%     original cuts until there is no more bias remaining and until at
%     least 9N more cuts have been made and the CV score stops changing by 
%     at least 10^-6 eV. 
%     This turns out to be a far more stable (repeatable) value than the
%     previous randomaly biased version.
%
%     2018-02-15: Optimized. Now allows for a 0.001% bias where bias is
%     defined as total imbalance divided by total number of cuts of the
%     data. Also skips a step if a matrix is 'close to singular'

if ~exist('optional_reference_CV_score','var')
    optional_reference_CV_score = 10;
end
ref_CV_score = optional_reference_CV_score;
lastwarn('')
warning('off','all')

% remove 0 coverage structure (the one with all zeros) from consideration


N = size(test_mat,1);
sz_constn_set = floor((1-a)*N);
sz_validn_set = ceil(a*N);


% calculate CV score for this CE

count=zeros(1,N);
diffCV = 1000;
new_CV_score = 100;
k = 1;
cnt_rng = 100;

tot_skips = 0;

while abs(diffCV) > 10^-7 || k < 10*N || cnt_rng/k > 0.0001
    old_CV_score = new_CV_score;

    if tot_skips > N
        CV_score = 10000;
        std_dev = 10000;
        return
    end
    
    if k <= N
        cc = randperm(N,sz_constn_set);
        constn_configs = cc;
        constn_set = test_mat(constn_configs,:);
        validn_configs = 1:N;
        validn_configs(constn_configs)=[];
        count(validn_configs)=count(validn_configs)+1;
        cnt_rng = max(count) - min(count);
    else
        [~, ii] = sort(count);
        validn_configs = ii(1:sz_validn_set);       
        constn_configs = 1:N;
        constn_configs(validn_configs)=[];
        constn_set = test_mat(constn_configs,:);
        count(validn_configs)=count(validn_configs)+1;
        cnt_rng = max(count) - min(count);     
    end
    constn_en = energy(constn_configs);
    validn_set = test_mat;
    validn_set(constn_configs,:) = [];
    validn_en = energy;
    validn_en(constn_configs) = [];
    [~,warcheck] = lastwarn;
    if strcmp(warcheck,'MATLAB:nearlySingularMatrix')
        lastwarn('');
        tot_skips = tot_skips + 1;
        continue
    end
    
    if any(all(constn_set==0,1))==1 
        tot_skips = tot_skips + 1;
        continue
    end
    coeffs = (constn_set'*constn_set)\(constn_set'*constn_en);
    pred_en = validn_set*coeffs;
    predn_err = pred_en - validn_en;
    sq_predn_err(k) = dot(predn_err,predn_err); %#ok<AGROW>
    sum_predn_err(k) = sum(predn_err); %#ok<AGROW>
    new_CV_score = sqrt(mean(sq_predn_err)/sz_validn_set);  
    diffCV = new_CV_score - old_CV_score;
    
    k = k + 1;
    if tot_skips > 1000
        new_CV_score = 10000;
        break
    end
    if new_CV_score - ref_CV_score > 0.0005  && abs(diffCV) < 0.001 && k > 2*N       
        break
    end
    if new_CV_score - ref_CV_score > 0  && abs(diffCV) < 0.01 && k > 10*N 
        break
    end
    if k > 500*N
        break
    end
    if new_CV_score > 1E6
        break
    end
end
%fprintf("...\n")
CV_score = new_CV_score;
avg_error = sum(sum_predn_err)/(k*sz_validn_set);
std_dev = sqrt(sum(sq_predn_err)/(k*sz_validn_set)-avg_error^2);


        

