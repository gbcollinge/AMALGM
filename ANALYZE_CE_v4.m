%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYZE_CE.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to hopefully find the optimal set of clusters that best
% describes and predicts the training data provided in matrix 
% "surf_normalized_counts" and either "surf_energy" or "mf_resids". 

% Must have CV_calc_v3.m in the active directory.

% v2: Amongst non-unique structures, makes sure the lowest energy duplicate
% is used in analysis
% v3: uses CV_calc_v3.m to determine errors based on variance in the ECIs
% over separate cuts of the data.
% v4: outputs confidence intervals using student t distribution
% (needed function from statistics toolbox: "tinv(p,nu)")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short g

%%%% Find all unique structures and amongst duplicates, select the
%%%% structure corresponding to the lowest energy inputted.
uCovs = unique(coverages);
sMat = surf_normalized_counts;
testing = unique(surf_normalized_counts,'rows','stable');
indMat = 1:size(surf_normalized_counts,1);
sEn = surf_energy;
totoss = [];
for ii = uCovs'
    sub_mat= sMat(abs(coverages-ii)<0.001,:);
    sub_en = sEn(abs(coverages-ii)<0.001);
    sub_ind = indMat(abs(coverages-ii)<0.001);
    [uSubMat,jja,jjb] = unique(sub_mat,'rows','stable');
    sz_dup = numel(jjb)-numel(jja);
    if sz_dup == 0
        continue
    end
    [count,~,idxcount] = histcounts(jjb,numel(jja));
    kk = 0;
    count(count>1) = 2;
    for jj = 1:numel(count)
        if count(jj)>1
            count(jj) = count(jj) +kk;
            kk = kk + 1;
        end
    end
    what = count(idxcount);
    mincount = min(what(what>1));
    maxcount = max(what);
    for jj = mincount:maxcount
        test = find(what==jj);
        [~, iik] = min(sub_en(test));
        totoss = [totoss sub_ind(test(test~=test(iik)))]; %#ok<AGROW>   
    end
 
    
end
tokeep = setdiff(indMat,totoss);    
MAT = surf_normalized_counts(tokeep,:);
N=size(MAT,1);                                                             % total number of known structures
M=size(MAT,2);                                                             % total number of clusters
mod_coverages = coverages(tokeep);                                         % grab the associated coverages

% User provided "natural" vectors of the metal surface. Be sure to use
% the vectors that create an obtuse angle. For FCC(111) this would be 
% [1 ; 0 ; 0] and [-1/2 ; sqrt(3)/2 ; 0] instead of [1 ; 0; 0] and 
% [1/2 ; sqrt(3)/2 ; 0]. Technically these are the a, b, and c vectors of  
% p(1 x 1) unit cell of your surface. If you intend to use the z
% coordinate position in defining adsorption site locations (only important 
% if you have numerous adsorption sites and they aren't all on the same
% plane.)...i.e. you can easily define them as having the same z position), 
% you should provide the c vector as it appears in your POSCAR divided by 
% the 1st NN distance. Otherwise, you can (and should) leave it as 
% [0 ; 0 ; 1]. Mind, whatever length a unit vector within this coordinate 
% system is will be the "natural unit" used from here on out.
% NOTE: This is the only place where cartesian coordinates should be 
% encountered!

ux = [1 ; 0 ; 0];
uy = [0 ; 1 ; 0];
uz = [0 ; 0 ; 1];

% If you know which clusters correspond to 1-body interactions in
% equilibrium with some chemical potential, you may want to subtract off
% the electronic (DFT) energy from that cluster. Add it here, if so, where
% the first entry is the cluster ID number and the second entry is the
% energy 

% mu_elec(1,:) = [1 -14.778323];

% User specified ECI value(s). First entry is the cluster ID #, second entry
% is the ECI value (in eV).

% ECI_val(1,:) = [1 -1.687979080];

% What color to use for each site type?

color_type{1} = 'black'; %Subsurface
color_type{2} = 'blue'; %Surface

% What size radius to use for each site type?

radius_type(1) = 0.303; %Subsurface
radius_type(2) = 0.303; %Surface

% If you'd like to use the mean field residuals as your fitting energies,
% please set "mf_flag" to "1". Else, set to "0". Default is "0"
                
mf_flag = 0;

% If you'd like to use the adsorption energy to fit the energies, set
% "ads_flag" to "1". Else, set to "0". Default is "0"

ads_flag = 1;

% User defined percentage of total number of known energies that become
% validation set during CV calculation (for a constant leave-x-out, set "a"
% to x/N...for leave-one-out, this would be 1/N)

a = 0.6;

% If you'd like to specify that certain clusters be added at the start (no
% guarantee the algorithm won't remove them, mind you) add them here. Note:
% if you don't want ANY starting clusters, just delete the numbers and
% leave an empty set...do not comment this out.

CE = [1     2     6    16   110   235   251   388   521   691   728];
    
% Do you want to check for unrepresented clusters? ("1" = yes, "0" = no)

check_flag = 0;

% Do you want to plot all the clusters found?
% "0" = no; "1" = yes; "2" = plot AND save to a png file

plot_clusters = 0;

% For the purposes of determining the minimum surfaces size for MC
% simulations, define a maximum NN distance that should count+-

maxNN = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deficient = (all(MAT(:,CE)==0,1));
deficient_clusters = CE(deficient);

if any(deficient)
    CE(deficient) = [];
    fprintf("There are no non-zero terms in the counts matrix for the following cluster(s)\n")
    disp(deficient_clusters)
    fprintf("This(these) clusters will be removed.\n")
    
end

% Grab appropraite energies
if mf_flag == 1
    energy = mf_resids;
else
    energy = surf_energy(tokeep);
end                       
    
if ads_flag == 1 && mf_flag == 0
    energy = ads_energy(tokeep);
    
    MAT = MAT./mod_coverages;
    MAT(~isfinite(MAT)) = 0;
else
    energy = surf_energy(tokeep);
end


% create random vector of interactions that will be used as a starting CE
% in addition to specified starting CE
out_CE = setdiff(1:M,CE);       % which clusters aren't part of the starting clusters?

CE = sort(CE); 
CE_orig = CE;
test_mat = MAT(:,CE);

if size(CE,2) >= rank(MAT)
    error('There are not enough known structures for this many starting clusters. Reduce the amount or add more structures and try again.')
end
fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\nThe starting Cluster Expansion (CE) is:\n')
disp(CE)
[CV_score,std_dev,ECI_errors] = CV_calc_v3(test_mat,energy,a);        %calulate the intial CV score


fprintf('Its CV score is: %8.6g eV/site\n',CV_score)
fprintf('Its standard deviation is: %8.6g eV/site\n',std_dev)
fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')

fID = fopen('CLUSTER_POSITIONS.txt');                           
tline = fgetl(fID);
line=0;
while ischar(tline)
    line = line + 1 ;
    tline = fgetl(fID);  
    if line == 1
        cols = cell2mat(textscan(tline, '%f'))';
    end
end
line = line - 1;

ncols = size(cols,2);
nMax= round(-7/2+1/2*sqrt(49+8*ncols));
nRs = nchoosek(nMax,2);
fclose(fID);
CLUSTER_POS = zeros(line,ncols);
INT_sz = line;
fID = fopen('CLUSTER_POSITIONS.txt');
tline = fgetl(fID);
for ii = 1:line
    tline = fgetl(fID);
    CLUSTER_POS(ii,:) = cell2mat(textscan(tline, '%f'));
end
fclose(fID);

interactions = INTERACTIONS(CE,:); 

final_mat = MAT(:,CE);
final_coeffs = (final_mat'*final_mat)\(final_mat'*energy);
final_en = final_mat*final_coeffs;
final_resids = energy - final_en;
resids_squared = final_resids.^2;
RMSE = sqrt(mean(resids_squared));
outlie = find(final_resids > 2*RMSE);
structures_to_check = atat_names(outlie);
final_coeffs_orig = final_coeffs;    


%v4: find confidence intervals
covar_mat = inv(final_mat'*final_mat);
variance = diag(covar_mat);
standard_error = sqrt(variance);
deg_free = size(final_mat,1) - size(final_mat,2);
t_dist = tinv(0.95,deg_free);
internal_error = sqrt(sum(final_resids.^2)/(deg_free));
confidence_intervals = t_dist.*standard_error*internal_error;

% aaa=0;
% 
% mod_mat = final_mat - ((2/(2+aaa)).*mean(final_mat));
% mod_coeffs = (mod_mat'*mod_mat)\mod_mat'*energy;
% mod_en = final_mat*mod_coeffs;
% mod_resids = energy - mod_en;
% 
% 
% ads_en = energy./coverages(iaa);
% ads_en(~isfinite(ads_en)) =[];
% mod_ads_en = mod_en./coverages(iaa);
% mod_ads_en(~isfinite(mod_ads_en)) =[];
% mod_cov = coverages(iaa);
% mod_cov(mod_cov == 0) = [];
% 
% scatter(mod_cov,ads_en);
% hold on
% scatter(mod_cov,ads_en,'r+')
% hold on 
% scatter(mod_cov,mod_ads_en,'b*')
% hold off
% drawnow

%scatter(coverages(iaa),final_resids,'r+')
%hold on
%scatter(coverages(iaa),mod_resids,'b*')
%hold off
%drawnow

% Find estimated ECI error
% tot_error = 1000;
chinew = 0;
% sqrd_mij = final_mat.^2;
% while abs(tot_error) > 0.0000001
%     chiold = chinew;
%     num_rep = sum(final_mat~=0,2);
%     clus_rep = sum(final_mat~=0,1);
%     stepA = final_resids./final_mat;
%     stepA2 = chiold*stepA;
%     stepA2(~isfinite(stepA2)) = 0;
%     stepB = stepA2.^2;
%     stepC = sum(stepB,1);
%     stepD = stepC./clus_rep;
%     stepE = stepD.*sqrd_mij;
%  
%     MSR_guess =sum(sum(stepE,2),1)/N;
% 
%     tot_error = RMSE^2 - MSR_guess;
%     frac_error = tot_error/RMSE^2;
%     
%     chinew = chiold + frac_error*sqrt(tot_error);
% end
 chi = chinew;
% intermA = (chi*final_resids./final_mat).^2;
% intermA(~isfinite(intermA)) = 0;
% intermB = sum(intermA,1)./clus_rep;
Estimated_ECI_errors = ECI_errors;
pcent_error = 100*Estimated_ECI_errors./abs(final_coeffs_orig');
structure_weight = chi;

% For plotting purposes
surf_x = ux;
surf_y = uy;
surf_x(3) = [];
surf_y(3) = [];
Us = [surf_x surf_y];
Cth = dot(surf_x,surf_y)/(norm(surf_x)*norm(surf_y));   % cos(theta) in terms of dot products
Sth = sqrt(1-Cth^2);
Rth = [Cth -Sth ; Sth Cth];                             % Rotation matrix

% Need the combinations of 2-body clusters for checking higher-body
% clusters
refcomb = [(1:nchoosek(nMax,2))' sortrows(nchoosek(1:nMax,2),2)];

if check_flag == 1
    % Check for un-represented clusters...ask the user how to deal with them
    if any(~any(surf_normalized_counts,1))          
        [row_INT, col_INT] = size(INTERACTIONS);                             % get the No. of rows and cols in the INTERACTION matrix
        nbody = round(1/2*(1+8*col_INT)^(1/2)-1/2);                          % find the original No. of n-bodies from the total number of columns
        unrep = find(~any(surf_normalized_counts,1));                        % find which columns (which correspond to clusters) in surf_normalized_counts are all zeros; these clusters aren't represented in the provided clusters
        sad_clusters = INTERACTIONS(unrep,:);                                % extract these from INTERACTIONS: call them "sad_clusters" (bc they are)
        unrepNN = unrep(sad_clusters(:,nbody+2)==0 & sad_clusters(:,2)~=0);  % find which of these clusters are of the 2-body type
        NN = INTERACTIONS(unrepNN,:);                                        % extract these from sad_clusters
        unrepElse = setdiff(unrep,unrepNN);                                  % find which clusters remain after extracting the 2-body types
        Else = INTERACTIONS(unrepElse,:);                                    % extract these from sad_clusters
        sad_R = NN(:,nbody+1);                                               % find the R distance  (R^2 in reality) of these 2-body clusters: "sad_R"
        fprintf('%d Clusters were found that are not represented in any of the provided structures.\n',size(unrep,2))
        fprintf('\nYou need to decide what to do with these clusters.\n\nType "inf" if this cluster is non allowed to occur during Monte Carlo.\nType "0" if this cluster is energy-neutral compared to the 0-coverage ads E (default).\nType an interaction energy (in units of eV) if you know what it should be (unlikely).\nOtherwise, you need to provide structures with this cluster in them!\n\n') 
        % loop through the 2-body clusters.
        dumCE = [];
        for j = 1:size(sad_R,1)                                          
            intrgate = NN(j,:);
            fprintf('How should this cluster be represented?\n')
            dumstring = sprintf('%6.4g\t', intrgate);
            fprintf('%s\n',dumstring)

            % plot the cluster

            the_types = CLUSTER_POS(unrepNN(j),1:nMax);
            nbodies = sum(the_types ~= 0);

            X = zeros(nbodies,1);
            Xdir = X;
            Y = X;
            Ydir = X;

            for jj = 1:nbodies
                Xdir(jj) = CLUSTER_POS(unrepNN(j),nMax+nRs+3*(jj-1)+1);
                Ydir(jj) = CLUSTER_POS(unrepNN(j),nMax+nRs+3*(jj-1)+2);
                Cart = Us*[Xdir(jj);Ydir(jj)];
                X(jj) = Cart(1);
                Y(jj) = Cart(2);
            end

            Xmin = floor(min(Xdir(:)))-1;
            Ymin = floor(min(Ydir(:)))-1;
            Xmax = ceil(max(Xdir(:)))+1;
            Ymax = ceil(max(Ydir(:)))+1;
            Xver = [Xmin:Xmax ; Xmin:Xmax];
            Yver = [Ymin*ones(1,(Xmax-Xmin+1)); Ymax*ones(1,(Xmax-Xmin+1))];
            Xhor = [Xmin*ones(1,(Ymax-Ymin+1)); Xmax*ones(1,(Ymax-Ymin+1))];
            Yhor = [Ymin:Ymax ; Ymin:Ymax];
            ver_bot = [Xver(1,:);Yver(1,:)];
            ver_top = [Xver(2,:);Yver(2,:)];
            hor_bot = [Xhor(1,:);Yhor(1,:)];
            hor_top = [Xhor(2,:);Yhor(2,:)];
            new_ver_bot = Us*ver_bot;
            new_ver_top = Us*ver_top;
            new_hor_bot = Us*hor_bot;
            new_hor_top = Us*hor_top;

            Xver = [new_ver_bot(1,:);new_ver_top(1,:)];
            Yver = [new_ver_bot(2,:);new_ver_top(2,:)];
            Xhor = [new_hor_bot(1,:);new_hor_top(1,:)];
            Yhor = [new_hor_bot(2,:);new_hor_top(2,:)];

            figure
            plot(Xver,Yver,'k')
            hold on 
            plot (Xhor,Yhor,'k')
            cluslab = ['Cluster ' num2str(unrepNN(j))];
            xlabel(cluslab);
            hold on
            for jj = 1:nbodies
                type = the_types(jj);
                circles(X(jj),Y(jj),radius_type(type),'facecolor',color_type{type});
                hold on
            end
            axis equal
            hold off
            drawnow

            % Ask what to do with them
            sad_E = input('"inf", "0", OR an interaction energy (in units of eV): ');  
            % ...and if that choice should be applied to subsequent clusters that are made up of this 2-body cluster
            alsoNbody = input('Apply this to all clusters that have this 2-body cluster in them?\n "yes" = "1" ; "no" = "0" (default): '); 
            if alsoNbody == 1                                   % yes, apply this change
                B1 = NN(j,1);                                       % body #1's type
                B2 = NN(j,2);                                       % body #2's type
                extraCE = [];
                for dd = 1:size(refcomb,1)
                    pass1 = Else(abs(Else(:,nMax+dd)-sad_R(j)) < 0.05,:);                                        % which clusters match the R  at R-position "dd"
                    pass2 = pass1((pass1(:,refcomb(dd,2)) == B1 & pass1(:,refcomb(dd,3)) == B2) | (pass1(:,refcomb(dd,2)) == B2 & pass1(:,refcomb(dd,3)) == B1),:);       % of these, which clusters have the correct body types                   
                    [~,newCE] = intersect(INTERACTIONS,pass2,'rows');                                           % find the correct clusters number(s) 
                    extraCE = [extraCE ; newCE]; %#ok<AGROW>                                                    % concatenate as we find the "extra"clusters

                end
            else                                                % no, don't apply this change
                extraCE =[];                                        % just make extraCE blank     
            end
            if sad_E ~= 0                                       % if sad_E IS "zero" we don't need to add this cluster--or any "extraCE"--to the CE
                [~,thisC] = intersect(INTERACTIONS,NN(j,:),'rows');
                extraCE = [ thisC extraCE' ];                      % concatentate this 2-body cluster and "extraCE" 
                extraEIC = ones(size(extraCE,2),1).*sad_E;          % create a repeated list of user's choice of EIC ("inf" or a non-zero interaction energy)

                % add the new clusters and EIC's to the CE and EICs ("final_coeffs)
                dumCE = [dumCE extraCE]; %#ok<AGROW>
                CE = [CE extraCE]; %#ok<AGROW>    

                final_coeffs = [final_coeffs ; extraEIC]; %#ok<AGROW>
            else
                [~,thisC] = intersect(INTERACTIONS,NN(j,:),'rows');
                dumCE = [dumCE thisC extraCE']; %#ok<AGROW>
            end
        end
        unrepElse = setdiff(unrep,dumCE); % Find which clusters remain that need addressing
        Else = INTERACTIONS(unrepElse,:); % Extract the clusters from INTERACTIONS
        % Now go through any remaining clusters that aren't NN ("Else") that
        % have survived the above process. Ask user if they'd like to apply a
        % single entry first
        fprintf('\nThere are %d clusters remaining',size(unrepElse,2))
        justdoit = input('\nThe remaining clusters are 3-body and above.\n...Would you like to apply a single value to all of them now? (enter "1")\n...or would you rather go through them one-by-one? (enter "0")\nYour choice:');
        skip_flag = 0;
        if justdoit == 1
            sad_E = input('"inf", "0", OR an interaction energy (in units of eV):\n');     
        end
        for j = 1:size(Else,1)
            intrgate = Else(j,:);
            if justdoit ~= 1
                fprintf('How should this cluster be represented?\n')
                dumstring = sprintf('%6.4g\t', intrgate);
                fprintf('%s\n',dumstring)

                % plot the cluster

                the_types = CLUSTER_POS(unrepElse(j),1:nMax);
                nbodies = sum(the_types ~= 0);

                X = zeros(nbodies,1);
                Xdir = X;
                Y = X;
                Ydir = X;

                for jj = 1:nbodies
                    Xdir(jj) = CLUSTER_POS(unrepElse(j),nMax+nRs+3*(jj-1)+1);
                    Ydir(jj) = CLUSTER_POS(unrepElse(j),nMax+nRs+3*(jj-1)+2);
                    Cart = Us*[Xdir(jj);Ydir(jj)];
                    X(jj) = Cart(1);
                    Y(jj) = Cart(2);
                end

                Xmin = floor(min(Xdir(:)))-1;
                Ymin = floor(min(Ydir(:)))-1;
                Xmax = ceil(max(Xdir(:)))+1;
                Ymax = ceil(max(Ydir(:)))+1;
                Xver = [Xmin:Xmax ; Xmin:Xmax];
                Yver = [Ymin*ones(1,(Xmax-Xmin+1)); Ymax*ones(1,(Xmax-Xmin+1))];
                Xhor = [Xmin*ones(1,(Ymax-Ymin+1)); Xmax*ones(1,(Ymax-Ymin+1))];
                Yhor = [Ymin:Ymax ; Ymin:Ymax];
                ver_bot = [Xver(1,:);Yver(1,:)];
                ver_top = [Xver(2,:);Yver(2,:)];
                hor_bot = [Xhor(1,:);Yhor(1,:)];
                hor_top = [Xhor(2,:);Yhor(2,:)];
                new_ver_bot = Us*ver_bot;
                new_ver_top = Us*ver_top;
                new_hor_bot = Us*hor_bot;
                new_hor_top = Us*hor_top;

                Xver = [new_ver_bot(1,:);new_ver_top(1,:)];
                Yver = [new_ver_bot(2,:);new_ver_top(2,:)];
                Xhor = [new_hor_bot(1,:);new_hor_top(1,:)];
                Yhor = [new_hor_bot(2,:);new_hor_top(2,:)];

                figure
                plot(Xver,Yver,'k')
                hold on 
                plot (Xhor,Yhor,'k')
                cluslab = ['Cluster ' num2str(unrepElse(j))];
                xlabel(cluslab);
                hold on
                for jj = 1:nbodies
                    type = the_types(jj);
                    circles(X(jj),Y(jj),radius_type(type),'facecolor',color_type{type});
                    hold on
                end
                axis equal
                hold off
                drawnow

                % Ask what to do with them
                sad_E = input('"inf", "0", OR an interaction energy (in units of eV):\n');           
            end
            if sad_E ~= 0                                   % if sad_E is "zero" we don't need to add this cluster to the CE          
                % add the new cluster and EIC to the CE and EICs ("final_coeffs)
                CE = [CE unrepElse(j)];         %#ok<AGROW>
                final_coeffs = [final_coeffs ; sad_E]; %#ok<AGROW>
            end
        end
    end
end

finalCE_pos = CLUSTER_POS(CE,:);


kk=1;
hh = 1;
q = ones(1,size(CE,2));
if plot_clusters == 1 || plot_clusters == 2
    flag1 = input('\nReady to plot the clusters in the CE.\nDo you want to plot the clusters with "inf" EIC?\n"no" = "0" ; "yes" = "1"\nChoice: ');
end
for ii = CE
    the_types = CLUSTER_POS(ii,1:nMax);
    nbodies = sum(the_types ~= 0);
    
    %if nbodies == 1
    %    continue
    %end
        
    X = zeros(nbodies,1);
    Xdir = X;
    Y = X;
    Ydir = X;
    
    for jj = 1:nbodies
        Xdir(jj) = CLUSTER_POS(ii,nMax+nRs+3*(jj-1)+1);
        Ydir(jj) = CLUSTER_POS(ii,nMax+nRs+3*(jj-1)+2);
        Cart = Us*[Xdir(jj);Ydir(jj)];
        X(jj) = Cart(1);
        Y(jj) = Cart(2);
    end
    
    if nbodies == 2 && length(unique(the_types(the_types~=0))) == 1 && ~isinf(final_coeffs(kk))               % If this is a 2-body (1st cond), homogeneous (2nd cond), with finite EIC (3rd condition) cluster
        ru = Xdir(2) - Xdir(1);
        rv = Ydir(2) - Ydir(1);
        rvec = [ru;rv];                         % created the vector (in "natural coordinates"), "rvec", of this 2-body cluster
        svec = Us\Rth*Us*rvec;                  % create the complimentary vector, "svec", that together with rvec produces the ordered structure correspoding to this 2-body cluster
        q(hh) = abs(det([rvec svec]));
        hh = hh + 1;
    end
    if plot_clusters == 1 || plot_clusters == 2
        if flag1 == 0
            if isinf(final_coeffs(kk))
                continue
            end
        end          
        Xmin = floor(min(Xdir(:)))-1;
        Ymin = floor(min(Ydir(:)))-1;
        Xmax = ceil(max(Xdir(:)))+1;
        Ymax = ceil(max(Ydir(:)))+1;
        Xver = [Xmin:Xmax ; Xmin:Xmax];
        Yver = [Ymin*ones(1,(Xmax-Xmin+1)); Ymax*ones(1,(Xmax-Xmin+1))];
        Xhor = [Xmin*ones(1,(Ymax-Ymin+1)); Xmax*ones(1,(Ymax-Ymin+1))];
        Yhor = [Ymin:Ymax ; Ymin:Ymax];
        ver_bot = [Xver(1,:);Yver(1,:)];
        ver_top = [Xver(2,:);Yver(2,:)];
        hor_bot = [Xhor(1,:);Yhor(1,:)];
        hor_top = [Xhor(2,:);Yhor(2,:)];
        new_ver_bot = Us*ver_bot;
        new_ver_top = Us*ver_top;
        new_hor_bot = Us*hor_bot;
        new_hor_top = Us*hor_top;

        Xver = [new_ver_bot(1,:);new_ver_top(1,:)];
        Yver = [new_ver_bot(2,:);new_ver_top(2,:)];
        Xhor = [new_hor_bot(1,:);new_hor_top(1,:)];
        Yhor = [new_hor_bot(2,:);new_hor_top(2,:)];

        figure
        plot(Xver,Yver,'k')
        hold on 
        plot (Xhor,Yhor,'k')
        cluslab = ['Cluster ' num2str(ii)];
        EIClab = ['EIC is ' num2str(final_coeffs(kk))];
        xlabel({cluslab , EIClab});
        hold on
        for jj = 1:nbodies
            type = the_types(jj);
            circles(X(jj),Y(jj),radius_type(type),'facecolor',color_type{type});
            hold on
        end
        axis equal
        hold off
        drawnow
        if plot_clusters == 2
            saveas(gcf, cluslab,'png')
        end
    end
    kk = kk + 1;
end

%%% Modify CE and EICs: Add back in specified cluster ID #s and EICs
if exist('ECI_val','var')==1
    if ~isempty(ECI_val)
        sz_ECI_vals = size(ECI_val,1);

        for ii = 1:sz_ECI_vals
            cluster_num = ECI_val(ii,1);
            cluster_ECI = ECI_val(ii,2);

            CE = [cluster_num CE]; %#ok<AGROW>
            energy = energy + orig_surf_normalized_counts(tokeep,cluster_num)*cluster_ECI;
            final_en = final_en + orig_surf_normalized_counts(tokeep,cluster_num)*cluster_ECI;
            final_coeffs = [cluster_ECI ; final_coeffs]; %#ok<AGROW>
        end
    end
end
[CE, iCE] = sort(CE);
final_coeffs = final_coeffs(iCE);

% Find the minimum MC surface size, write to file
q(abs(round(q)-q) > 0.0001) = [];      % remove non-integer values
qwmaxNN = q(q<=maxNN);
if qwmaxNN > 1
    Nsx_min = lcms(qwmaxNN);
else
    Nsx_min = 1;
end

fID = fopen('Nsx_min.txt','wt+');
fprintf(fID,'NN distances for each 2-body cluster:\n');
fclose(fID);
dlmwrite('Nsx_min.txt',q(abs(q-1)>0.0001),'delimiter',' ','-append')
fID = fopen('Nsx_min.txt','A');
fprintf(fID,'You have specified a maximum NN distance of\n%d\n', maxNN);
fprintf(fID,'So the minimum surface length is:\n');
fclose(fID);
dlmwrite('Nsx_min.txt',Nsx_min,'-append')

EICs_w_clusnum = [CE_orig' final_coeffs_orig Estimated_ECI_errors' pcent_error' confidence_intervals];

mean_ECI_errors = mean(Estimated_ECI_errors);
mean_confidence_intervals = mean(confidence_intervals);


% Find the minimum sublattice size
minSL = ceil(max(q(:)));
dlmwrite('min_sublat.txt',minSL)

%%% Write CE to MC_POSITIONS folder
if exist('MC_POSITIONS','dir')==7
    cd './MC_POSITIONS'
    dlmwrite('CE',CE')
    cd '../'
end

%%% Write ECIs to MC_POSITIONS folder
if exist('MC_POSITIONS','dir')==7
    cd './MC_POSITIONS'
    dlmwrite('ECIs',final_coeffs)
    cd '../'
end

%%% Write final result to file
dlmwrite('LG_Model.txt',EICs_w_clusnum,'delimiter','\t')

%%% Write the output to file "LG_EICs.txt"
dlmwrite('LG_EICs.txt',CE,'delimiter','\t')
dlmwrite('LG_EICs.txt',final_coeffs,'delimiter','\t','-append')
fID = fopen('LG_EICs.txt','a');
fprintf(fID,'%8.6f eV/site',CV_score);
fclose(fID);

%%% Write coverages and surface energies and predicted surface energies
all_data = [mod_coverages energy final_en final_resids];
dlmwrite('all_data.txt',all_data,'delimiter','\t')

%%% User Displayed Outputs
fprintf('The LG ECI for this CE are:\n')
fprintf('Cluster#\tECI     \tECI_error\t  percent_error confidence intervals (eV)\n')
disp(EICs_w_clusnum)

disp([mean_ECI_errors mean_confidence_intervals mean_ECI_errors-mean_confidence_intervals])

fprintf('The RMSR is:\n')
disp(RMSE)
fprintf('Structures with residuals in excess of 2*RMSR:\n')
disp(structures_to_check)
fprintf('This CE corresponds to the following interactions:\n')
disp(interactions)
fprintf('Type "final_resids" to see what the corresponding residuals are:\n')

fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
clearvars -except CE slope_counts slope_energy part_slope_en true_gs_slopes_counts true_gs_slopes part_ads_en ads_normalized_counts true_gs iso_gs diff_ads_counts diff_surf_counts diff_surf_energies diff_ads_energies true_gs_covs true_gs_ads_e  true_gs_ens true_gs_slopes ads_energy part_en orig_surf_normalized_counts atat_names newest_configs residuals mf_resids counts surf_normalized_counts coverages surf_energy sMax INTERACTIONS

