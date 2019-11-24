0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CE_Generator.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build and find interaction terms in a cluster expansion for adsorbates
% on metal surfaces, along with occupancies of these interactions (given a
% certain configuration)

% v5:   finds "problem structures"
% v6:   find which of the new structures in "new_configs" is repeated, if
% any; does not include these in the calculation of the external CV score.
% v6.1: fixes "double counting" of interactions based on site number
% v7: finds ground states
% v8: determines the matrix of constraints needed to ensure that
% the ground states are recovered in the fitting of the ECIs to the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('surf_normalized_counts','var') == 1
    MAT = surf_normalized_counts;
    EN = surf_energy;
end
clearvars -except MAT EN surf_normalized_counts surf_energy
format short g

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User provided "natural" vectors of the metal surface. Be sure to use
% the vectors that create an obtuse angle. For FCC(111) this would be 
% [1 ; 0 ; 0] and [-1/2 ; sqrt(3)/2 ; 0] instead of [1 ; 0; 0] and 
% [1/2 ; sqrt(3)/2 ; 0]. Technically these are the a, b, and c vectors of  
% p(1 x 1) unit cell of your surface. If you intend to use the z
% coordinate position in defining adsorption site locations (only important 
% if you have numerous adsorption sites and they aren't all on the same
% plane.)...i.e. you can easily define them as having the same z position), 
% you should provide the c vector as it appears in your POSCAR. Otherwise, 
% you can (and should) leave it as [0 ; 0 ; 1]. Mind, whatever length a
% unit vector within this coordinate system is will be the "natural unit" 
% used from here on out.
% NOTE: This is the only place where cartesian coordinates should be 
% encountered!

ux = [1 ; 0 ; 0];
uy = [-1/2 ; sqrt(3)/2 ; 0];
uz = [0 ; 0 ; 1];

% Alternatively, change infile to "1" and provide a file called
% "NATURAL_COORDINATES.txt" with each ux uy and uz provided as column
% vectors. 

infile = 0;

% This file must be written in decimal (floating point) format. Be careful
% here, the script appears to suffer from round off errors and you might
% need to a "fudge factor" to your Rmax. Experiment. Otherwise, don't use
% this feature.
% An example for FCC(11) or HCP(0001) follows:
%
%     1 0 0
%     -0.5 0.866025403784439 0
%     0 0 1    
%
%     end of example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you know which clusters correspond to 1-body interactions in
% equilibrium with some chemical potential, you may want to subtract off
% the electronic (DFT) energy from that cluster. Add it here, if so, where
% the first entry is the cluster ID number and the second entry is the
% energy 

mu_elec(1,:) = [1 -14.778323];

% User specified ECI value(s). First entry is the cluster ID #, second entry
% is the ECI value (in eV).

%ECI_val(1,:) = [1 -1.687979080];

% Specify whether to make regression constraints based on surface energy
% ("ads_flag = 0") or adsorption energy ("ads_flag = 1")

ads_flag = 1;

% How many sites within this "natural" cell?

sites_per_cell = 1;

% If there are more than one site, please provide their location within the
% natural unit cell (in "natural coordinates") along with a number to
% signify the type of site. The site "type" can be any positive non-zero
% integer and does not need to be continuous. 
% e.g. if you have an FCC(111) surface, there are potentially top sites (1),
% bridge sites (2), fcc hollow sites (3), and hcp hollow sites (4). If
% you want to specify more than one type of adsorbate, you can do that here
% by repeating the same adsorption site, but changing the "type" number 
% (i.e. the 4th element)
    Site(:,1) = [0; 0 ; 0 ; 1]; % Subsurface
   %Site(:,2) = [1/2; 1/2 ; 0 ; 2]; % Surface

% If any of the sites are linked, as in through a bond, then identify
% below. THis will simply remove the point EIC (V naught) for the linked
% site

linked = [];

% Which sites (not types) will be used to calculate the coverage of this
% system?

coverage_sites = [1];

% User specified overall maximum N-body clusters to include (even if the
% max is different for different site types, still specify the max of all
% types here)

maxNbody= 5;

% User specified "problematic length". When a supercell has a length that 
% is equal to or smaller than this, the corresponding structure will be 
% marked as problematic and REMOVED.

prob_length = 1;

%%%%%%%%%%%%%%%%%%% BOOK KEEPING, PLEASE DON'T TOUCH %%%%%%%%%%%%%%%%%%%%%%
        vecbody = [ones(1,maxNbody)*maxNbody maxNbody];
        Rmax = zeros(vecbody);
%%%%%%%%%%%%%%%%%%%%%%%%%% OKAY DONE, CONTINUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%      

% User defined maximum interaction distances "Rmax" in units of 
% natural unit vectors. Each matrix (e.g. Rmax(:,:,2)) corresponds to a 
% n-body interaction (e.g. 2 body interaction). Each row and column 
% correspond to each site type (so if there are 3 site TYPES, these
% will be 3 x 3 symmetric matrices. Each element corresponds to 
% interactions between the types designated by the row and column. For
% example, if there are 3 site types (1, 2, and 3), Rmax(:,:,3) contains
% the maximum site distances (or "cluster sizes") for 3-body interactions 
% and Rmax(1,3,3) is the maximum 3-body site distance between sites 
% 1 and 3 (corresponding to the sites entered above). If you want (say) 4
% body interactions between site 1 and itself (Rmax(1,1,4)) but not between
% site 1 and 2, just enter "0" for that entry (i.e. Rmax(1,2,4) = 0).

Rmax(1,1,:,:,:,2) = 5;
Rmax(1,1,1,:,:,3) = 3.7;
Rmax(1,1,1,1,:,4) = 3.7;
Rmax(1,1,1,1,1,5) = 3;

% User defined minimum interaction distance "Rmin" in units of 
% natural unit vectors. This is the same for all types of interactions.

Rmin = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%& END USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
display_flag = 0; 
% manipulate/fill in Rmax array
for ii=2:maxNbody
    Ncombs = nmultichoosek(unique(Site(4,:)),ii);
    sz_Ncombs = size(Ncombs,1);
    RtoUse = zeros(sz_Ncombs,1);
    for jj = 1:sz_Ncombs
        Nperms = unique(perms(Ncombs(jj,:)),'rows');
        sz_Nperms = size(Nperms,1);
        for kk = 1:sz_Nperms
            vecp = [Nperms(kk,:) ones(1,maxNbody-ii) ii];
            vecind = num2cell(vecp); 
            Rvecind(jj,kk) = Rmax(sub2ind(size(Rmax),vecind{:})); %#ok<SAGROW>
        end
        RtoUse(jj) = max(Rvecind(jj,:));
    end
    if any(RtoUse == 0)
        whichjj = find(RtoUse == 0);
        uniqtype = unique(Ncombs(whichjj,:));
        for jj = 1:sz_Ncombs
            if all(unique(Ncombs(jj,:)) == uniqtype)
                matchjj = jj;
                break
            end
        end
        RtoUse(whichjj) = RtoUse(matchjj);
    end
    for jj = 1:sz_Ncombs
        Nperms = unique(perms(Ncombs(jj,:)),'rows');
        sz_Nperms = size(Nperms,1);
        for kk = 1:sz_Nperms
            vecp = [Nperms(kk,:) ones(1,maxNbody-ii) ii];
            vecind = num2cell(vecp); 
            Rmax(sub2ind(size(Rmax),vecind{:})) = RtoUse(jj);
        end    
    end
end
Rmax = Rmax.^2;
Rmin = Rmin.^2;

fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
fprintf('Running COUNTS_GEN.m.\nAll configurations available should be placed in a folder name "configs".\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\nWorking...\n\n')
fprintf('----------------------------------------------------------------\n\n')

% Get the natural coordinate system from file "NATURAL_COORDINATES.txt" if
% infile flag has been turned on
if infile == 1
    fID = fopen('NATURAL_COORDINATES.txt');
    tline = fgetl(fID);
    ux = cell2mat(textscan(tline, '%f'));
    tline = fgetl(fID);
    uy = cell2mat(textscan(tline, '%f'));
    tline = fgetl(fID);
    uz = cell2mat(textscan(tline, '%f'));   
    fclose(fID);
end

% Determine the "norm conserving matrix" for later determination of
% distances
natcoor = [ux uy uz];
normR = natcoor'*natcoor;
    
% determine the max X and Y values needed to reach the maximum Rmax value
% specified
unitX = [1;0;0];
unitY = [0;1;0];
lengthX = unitX'*normR*unitX;
lengthY = unitY'*normR*unitY;
maxRmax = max(Rmax(:));
factor = lengthX/lengthY;
if factor > 1
    big_vec = factor*unitY + unitX;
    
    leng_bigvec = big_vec'*normR*big_vec;
    max_fac = sqrt(maxRmax)/leng_bigvec;
    
    maxX = ceil(max_fac*big_vec(1));
    maxY = ceil(factor*max_fac*big_vec(2));
elseif factor < 1
    big_vec = unitY + unitX./factor;
    
    leng_bigvec = big_vec'*normR*big_vec;
    max_fac = sqrt(maxRmax)/leng_bigvec;
    
    maxX = ceil(max_fac*big_vec(1)/factor);
    maxY = ceil(max_fac*big_vec(2));
else
    big_vec = unitY + unitX;
    
    leng_bigvec = big_vec'*normR*big_vec;
    max_fac = sqrt(maxRmax)/leng_bigvec;
    
    maxX = ceil(max_fac*big_vec(1));
    maxY = ceil(max_fac*big_vec(2));
end
maxX= ceil(maxX);
maxY=ceil(maxY);
% Total number of sites and bodies
Site = Site(:,Site(4,:)~=0);
sMax = size(unique(Site(4,:)),2);
nMax = maxNbody;
site_max = size(Site,2);

% Total number of body-to-body pair distances
numRs = nchoosek(nMax,2);

% Get all possible interactions from file "INTERACTIONS.txt" which 
% needs to be in the parent directory

fID = fopen('INTERACTIONS.txt');                           
tline = fgetl(fID);
line=0;
while ischar(tline)
    line = line + 1 ;
    tline = fgetl(fID);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Kill program if number of columns in INTERACTIONS.txt is inconsistent
    %%% with the number in those specified by nMax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if line == 2
        COLS = cell2mat(textscan(tline,'%f'))';
        INT_COLS = size(COLS,2);
        nCalc = -1/2 + 1/2*sqrt(1+8*INT_COLS);

        if abs(nCalc-nMax)>0.1
            error('Error. The maximum number of n-body interactions suggested by the number of columns in the INTERACTIONS.txt file is %3.0g, but you have specified %3.0g in this script. Correct this and try again.',nCalc,nMax)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
line = line - 1;

fclose(fID);
INTERACTIONS = zeros(line,nMax+numRs);
INT_sz = line;
fID = fopen('INTERACTIONS.txt');
tline = fgetl(fID);
for ii = 1:line
    tline = fgetl(fID);
    INTERACTIONS(ii,:) = cell2mat(textscan(tline, '%f'));
end
fclose(fID);

% Check if folder new_configs exists. If so, check that user has generated
% file LG_EICs.txt and/or mf_coeffs.txt. Prompt for user input where
% necessary.
if exist('new_configs','dir') == 7
    fprintf('Directory "new_configs" has been created...\n...checking if a LG model has been provided...\n') 
    if exist('LG_EICs.txt','file') == 0
        no_LG = 1;
        fprintf('"LG_EICs.txt" not found.\nWould you like this script to simply move any new structures\nfound in folder "new_configs" to folder "configs"?\n')
        move_input = input('\nEnter "1" (yes) or "0" (no): ');            
    else
        no_LG = 0;
        fprintf('"LG_EICs.txt" found.\nThis script will calculate the external-CV score for any new configurations\nfound in folder "new_configs".\n')
        fprintf('Do you want to move any non-problematic structures found in new_configs to folder "configs"?\n')
        move_input2 = input('\nEnter "1" (yes) or "0" (no): ');    
        fprintf('Do you want to move any problematic structures found in new_configs to folder "prob_configs"?\n')
        move_input3 = input('\nEnter "1" (yes) or "0" (no): ');    
    end
    if exist('mf_coeffs.txt','file') == 0
        no_mf = 1;
    else
        no_mf = 0;
        fprintf('\nmf_coeffs.txt found.\nThis script assumes the LG model in LG_EICs.txt is\nfit to the residuals of this mean field model.\nIf this is not the case,the external-CV score will be useless.\n')
        kill_flag = input('Proceed? ("1") or end this script? ("0"): ');
        if kill_flag ~= 1
            return
        end
    end
    fprintf('Proceeding...\n\n')   

    % Loop over configurations in directory "new configs", which should be located
    % in your working directory where this script is located
    all_configs = dir('./new_configs');             % command 'dir' assigns each file a fileID starting from 3
    newest_configs = {all_configs.name};
    newest_configs = newest_configs(3:end);
    num_configs = numel(all_configs)-2;             % since fileID 1 and 2 are STDIN and STDERR we want to discount these
    counts = zeros(num_configs,line);               % A whole bunch of preallocation follows...
    energy = zeros(num_configs,1);
    coverages = zeros(num_configs,1);
    surf_normalized_counts = counts;
    surf_energy = energy;
    confnum = 0;
    for jj=1:num_configs                    % start working through the new configurations...   
        if exist('move_input','var')== 1
            if move_input == 1 
                movefile(['./new_configs/' all_configs(jj+2).name],'./configs/')
                fprintf('NEW CONFIGURATION #%d, %s:\n',jj, all_configs(jj+2).name)
                if move_input2 ~= 1
                    fprintf('\tMoved!');              
                end
                continue
            end
            if move_input ~= 1
                fprintf('Ignoring new_configs.\n')
                break
            end
        end
        fprintf('NEW CONFIGURATION #%d, %s:\n',jj, all_configs(jj+2).name)
        backhome = cd('./new_configs/');
        confnum = confnum + 1;
        analyzed_flag = 0;
        
        % Determine how many lines there are in this file
        fID = fopen(all_configs(jj+2).name);
        tline = fgetl(fID);
        line=0;

        while ischar(tline)
            if isempty(tline)
                tline = fgetl(fID);
                continue
            end

            if strcmp(tline,'analyzed')                         % check if this configuration has already been analyzed
                analyzed_flag = 1;                              % Turn the analyzed_flag on
                tline = fgetl(fID);
                if strcmp(tline,'Confidence') || strcmp(tline,'No Confidence')
                    tline = fgetl(fID); %#ok<NASGU>
                end
                tline = fgetl(fID);                             % advance to where the counts should be               
                counts(jj,:) = cell2mat(textscan(tline, '%f')); % grab these counts for the matrix
                break
            else
            line = line + 1 ;
            tline = fgetl(fID);   
            end
        end   
        fclose(fID);    

        line = line -1;
        CELL = zeros(2,2);
        SIG = zeros(line-2,3);

        % Now extract the unit cell of the configuration (CELL) along with the
        % positions of the adsorbates (SIG) and the configuration energy
        % (energy)
        fID = fopen(all_configs(jj+2).name);

        for linenum = 1:2
            tline = fgetl(fID);
            CELL(linenum,:) = cell2mat(textscan(tline, '%f'));
        end
        for linenum = 3:line
            tline = fgetl(fID);
            SIG(linenum-2,:) = cell2mat(textscan(tline, '%f'));
        end
        tline = fgetl(fID);
        energy(jj) = str2double(tline);         

        fclose(fID);
        cd('../')


        Ns = norm(det([CELL(1,:);CELL(2,:)]))*sites_per_cell;
        rel_SIG = [];
        for ii = 1:max(size(coverage_sites))
            rel_SIG = [rel_SIG; SIG(SIG(:,3)==coverage_sites(ii),:)];  %#ok<AGROW>
        end
        num_adsorbates = size(rel_SIG,1);
        
        coverages(jj) = num_adsorbates/Ns;

        num_adsorbates = size(SIG,1);

        % find length of all effective vector lengths in the unit cell
        cell_length = zeros(1,4);
        cell_length(1) = CELL(1,:)*normR(1:2,1:2)*CELL(1,:)';
        cell_length(2) = CELL(2,:)*normR(1:2,1:2)*CELL(2,:)';
        cell_length(3) = (CELL(1,:)-CELL(2,:))*normR(1:2,1:2)*(CELL(1,:)-CELL(2,:))';
        cell_length(4) = (CELL(1,:)+CELL(2,:))*normR(1:2,1:2)*(CELL(1,:)+CELL(2,:))';
        
        % If any of those length are equal to or less than the problematic
        % length, flag this structure
        prob_flag = 0;
        if sum(cell_length-prob_length^2 < 1) > 0.999 && sum(cell_length-prob_length^2 < 1) < 1.999
            if exist('prob_configs','dir') ~= 7
                mkdir('prob_configs');
            end
            prob_flag = 1;            
        end
        
        if analyzed_flag == 1                                   % if this configuration has already been analyzed 
            fprintf('Already analyzed!\n')
            surf_normalized_counts(jj,:) = counts(jj,:)./Ns;
            surf_energy(jj) = energy(jj)./Ns;
            if move_input2 == 1
                movefile(['./new_configs/' all_configs(jj+2).name],'./configs/')
            end
            continue                                            % ...then we can skip the rest of this iteration
        end

        % Initial Size of sig matrix (dynamic preallocation of memory)
        BLOCK = 100;
        list_size = BLOCK;
        col_BLOCK = 4;
        sig = zeros(BLOCK,col_BLOCK);

        % Find SigMaxX and SigMaxY: the dimensions needed to create a surface
        % large enough to encompass the maxR distance for each adsorbate

        diagonal_vec = (CELL(1,:)+CELL(2,:));
        SigMaxX = diagonal_vec(1)+maxX;
        SigMaxY = diagonal_vec(2)+maxY;

        cellvecX = CELL'\[2*SigMaxX ; -2*SigMaxY];
        cellvecY = CELL'\[-2*SigMaxX ; 2*SigMaxY];

        xMin = cellvecY(1);
        xMax = cellvecX(1);

        if xMax < xMin
            dumX = xMax;
            xMax = xMin;
            xMin = dumX;
        end
        
        xMin = floor(xMin);
        xMax = ceil(xMax);

        yMin = cellvecX(2);
        yMax = cellvecY(2);

        if yMax < yMin
            dumY = yMax;
            yMax = yMin;
            yMin = dumX;
        end
        
        yMin = floor(yMin);
        yMax = ceil(yMax);

        xTot = xMax - xMin;
        yTot = yMax - yMin;

       % Populate "sig" matrix
        kk = 1;
        SIG_dif = zeros(num_adsorbates,2);
        for ii = 1:num_adsorbates     
            for xx = 0:xTot
                sigx = xx + xMin;         
                for yy = 0:yTot    
                    sigy = yy + yMin;
                    repeat = SIG(ii,1:2)+(CELL'*[sigx;sigy])';
                    SIG_dif = SIG(:,1:2) - repeat;                      % difference between this repeated position and all adsorbates' positions                   
                    R_check = (diag([SIG_dif zeros(num_adsorbates,1)]*normR*[SIG_dif zeros(num_adsorbates,1)]'));               
                    if any(R_check <= maxRmax+1)
                        Shift = Site(1:3,SIG(ii,3))';
                        sig_type = Site(4,SIG(ii,3));
                        sig(kk,:) = [[repeat 0]+Shift sig_type];
                        kk = kk + 1;                      
                    end
                    % Add new block of memory to sig matrix if needed
                    check = sig(any(sig~=0,2),:);
                    if size(check,1)/size(sig,1) > 0.95                 % Only 5% of the current allocation is left
                        list_size = size(check,1) + BLOCK;
                        sig(kk+1:list_size,:)=0;                        % ...so add a new block of memory
                    end                         
                end
            end  
        end

        sig = sig(any(sig~=0,2),:);
        sig = sortrows(sig,1:3);

        nTot = size(sig,1);
        X = sig(:,1);
        Y = sig(:,2);
        Z = sig(:,3);
        % Find where each site within the unit cell is within the sig matrix 
        sigPos = zeros(1,num_adsorbates);
        for ii = 1:num_adsorbates
            Shift = Site(1:3,SIG(ii,3))';
            sig_type = Site(4,SIG(ii,3));
            check = [[SIG(ii,1:2) 0]+Shift sig_type];
            [~,sigPos(ii)] = intersect(sig,check,'rows');
        end       

        % All combinations of pairs for every adsorbate on the SigMaxX by
        % SigMaxY surface
        CombTot = sortrows(nchoosek(1:nTot,2),2);
        nCombs = size(CombTot,1);
        frstn = CombTot(:,1)';
        scndn = CombTot(:,2)';

        % Find the delta X,Y,Z between each pair of bodies with each site
        delX = zeros(1,nCombs);
        delY = zeros(1,nCombs);
        delZ = zeros(1,nCombs);
        delXYZ = zeros(3,nCombs);
        R = zeros(1,nCombs);
        for ii = 1:nCombs
            delX(ii) = X(scndn(ii))-X(frstn(ii));
            delY(ii) = Y(scndn(ii))-Y(frstn(ii));
            delZ(ii) = Z(scndn(ii))-Z(frstn(ii));
            delXYZ(:,ii) = [delX(ii);delY(ii);delZ(ii)];
            R2 = delXYZ(:,ii)'*normR*delXYZ(:,ii);
            R(ii) = (R2);
        end

        % Associate R's with each of its associated pair's positions.
        % Swaping these positions changes nothing so do that now.
        Rassoc = zeros(nTot);
        for h = 1:size(CombTot,1)
            Rassoc(CombTot(h,1),CombTot(h,2)) = R(h);
            Rassoc(CombTot(h,2),CombTot(h,1)) = Rassoc(CombTot(h,1),CombTot(h,2));
        end

        % Now find all pairs that contain the number of the position of each
        % supercell adsorbate. Use this to extract all the R's that connect the
        % various n-body interactions    
        clusters = zeros(BLOCK, nMax+numRs);
        min_clusters = zeros(1,BLOCK);
        Int_sz = BLOCK;
        int_ptr = 1;

        % Create point (V naught) clusters first
        principle_sites = unique(Site(4,:));
        if size(principle_sites,2) > 1
            uu = unique(linked);
            nn = histcounts(linked);
            over_linked = uu(nn>1);
            principle_site = intersect(principle_sites,uu);

            for ii = 1:size(over_linked,2)
                over_linked_sites = linked(linked==over_linked(ii),:);
                principle_sites = [principle_sites min(over_linked_sites(:))]; %#ok<AGROW>
            end
            non_linked_sites = principle_sites;
            for ii = 1:size(linked,1)
                if isempty(intersect(linked(ii,:),principle_sites))
                    non_linked_sites = [non_linked_sites min(linked(ii,:))]; %#ok<AGROW>
                end
            end
        end

            non_linked_sites = principle_sites;
        for ii = 1:size(non_linked_sites,2)
       num_sites_here = size(SIG(SIG(:,3)==ii,:),1);
            for  ww= 1:num_sites_here
                clusters(int_ptr,1) = non_linked_sites(ii);
                posclusters(int_ptr,1) = non_linked_sites(ii);
                int_ptr = int_ptr + 1;
            end
        end

        % Now go through and find the rest of the clusters
        for ii = 2:nMax
            fprintf('\tFinding all %d-body clusters...\n',ii)
            for kk = 1:num_adsorbates
                fprintf('\t\tSite #%d\n',kk)
                curr_ads_type = sig(sigPos(kk),4);
                vec_dif = sig(:,1:3) - sig(sigPos(kk),1:3);
                R_vec_dif = (diag(vec_dif*normR*vec_dif'));       
                sigwR = [sig R_vec_dif];
                irrel_sigA = sigwR(sig(:,1)==sig(sigPos(kk),1) & sig(:,2)<sig(sigPos(kk),2),:);       % Find the "irrelevant" positions in sig (A) those where Y < X when X = the X of adsorbate kk
                irrel_sigB = sigwR(sig(:,1)<sig(sigPos(kk),1),:);                                     % ...(B) those where X is lower than the adsorbate        
                irrel_sigD = sigwR(sig(:,3)<sig(sigPos(kk),3),:);                                     % ...(D) whose site number is less than the adsorbates
            ind_sigC = zeros(10000,1);
            ind_ptr = 1; 
            for bb = 1:sMax                                                            % ...(C) those where their distance from the adsorbate is greater than this type's n-body Rmax                           
                nd = ndims(Rmax);                                                       % Since it can change, how man dimension in the Rmax array?
                hold_ind = repmat({':'},1,nd-3);                                        % Create nd-3 ":" indices to insert into the Rmax indexing
                subsetR = Rmax(curr_ads_type,bb,hold_ind{:},ii);                        % Get the 2-body subset of Rmax of the current site type and bb
                maxSubsetR = max(subsetR(:));                
                temp_ind = find(sigwR(:,4)==bb & sigwR(:,5) > maxSubsetR); 
                sz_temp = size(temp_ind,1);
                ind_sigC(ind_ptr:(ind_ptr+sz_temp-1)) = temp_ind;
                ind_ptr = ind_ptr + sz_temp;
            end          
            ind_sigC = ind_sigC(ind_sigC ~= 0,:);
            ind_sigC = sort(ind_sigC);
            irrel_sigC = sigwR(ind_sigC,:);                                                                                                               
            irrel_sig = [irrel_sigA ; irrel_sigB ; irrel_sigC; irrel_sigD];
                irrel_sig = setdiff(irrel_sig,sigwR(sigPos(kk),:),'rows');
                [rel_sig,rel_bodies] = setdiff(sigwR,irrel_sig,'rows');                               % Find which bodies (indices in sig) fit the above criteria
                if size(rel_bodies,1) <= 1
                    continue
                end
                rel_combs = nchoosek(rel_bodies',ii);                                               % Use these bodies to find the relevant n-body combos
                rel_combs = rel_combs(any(rel_combs == sigPos(kk),2),:);                            % Restrict the relevant bodies to those with the kk adsorbate in it
                bTot = size(rel_combs,1);
                dist = zeros(bTot,numRs);
                types = reshape(sig(rel_combs',4),ii,[])';                                          % rel_bodies are read row THEN column by sig, so transpose rel_bodies first, then read their types (col 4 of sig). Reshape this back into a matrix the same shape as rel_bodies
                good_types = zeros(bTot,ii);
                zz = 1;
                for ll = 1:bTot
                    pairs = sortrows(nchoosek(rel_combs(ll,:),2),2);
                    num_pairs = size(pairs,1);   
                    dist_temp = diag(Rassoc(pairs(:,1),pairs(:,2)))';
                    index_vec = num2cell([types(ll,:) ones(1,nMax-ii) ii]);
                    if all(dist_temp <= Rmax(sub2ind(size(Rmax),index_vec{:})))   % will need to generalize this part another time
                        dist(ll,1:num_pairs) = dist_temp; 
                        good_types(ll,:) = types(ll,:);  
                    end
                end                   
                fill_zeros = zeros(bTot,nMax-ii);
                dum_types = [good_types fill_zeros];              
                dum = [dum_types dist];  
                dum(all(dum == 0,2),:) = [];
                dum = dum(all(dum(:,nMax+1:nMax+num_pairs) >= Rmin-0.01,2),:);
                dum_sz = size(dum,1);                 
                    % Fill memory as needed
                    check = clusters(all(clusters == 0,2),:);
                    check_sz = size(check,1);               
                    while check_sz < dum_sz
                        Int_sz = Int_sz + BLOCK;
                        clusters(int_ptr+1:Int_sz,:)=0;
                        min_clusters(int_ptr+1:Int_sz) =0;
                        check = clusters(all(clusters == 0,2),:);
                        check_sz = size(check,1);               
                    end     

                clusters(int_ptr:int_ptr+dum_sz-1,:) = dum;
                int_ptr = int_ptr + dum_sz + 1;                     
            end 

        end 
        % Clean up clusters matrix: remove unused space
        clusters = clusters(any(clusters~=0,2),:);
        TotInts = size(clusters,1);

        % Find where each n-body interaction starts and ends
        section_ctr = zeros(1,nMax);
        dummy_int = clusters;
        for mm = 2:nMax
            if mm ~= nMax
                section = dummy_int(all(dummy_int(section_ctr(mm-1)+1:end,mm+1)==0,2),:);
            else
                section = dummy_int(section_ctr(mm-1)+1:end,:);
            end 
            section_ctr(mm) = section_ctr(mm-1) + size(section,1);
        end
        section_ctr(nMax+1) = TotInts;
        fprintf('\nAll n-body clusters found!\n\nNow assigning them to the provided interactions...\n')    

        % Now permute the n-bodies and check against the INTERACTIONS matrix

        for qq = 2:nMax
            Nperm = sortrows(perms(1:qq));        % Find all permutations of qq bodies
            Nperm_sz = size(Nperm,1);
            tpair = sortrows(nchoosek(1:qq,2),2);              % Total combinations of pairs used here according to number qq
            tpair_sz = size(tpair,1);
            for ii = section_ctr(qq-1)+1:section_ctr(qq)
                newRassoc = zeros(nMax);
                for kk = 1:tpair_sz
                    newRassoc(tpair(kk,1),tpair(kk,2)) = clusters(ii,nMax+kk);
                    newRassoc(tpair(kk,2),tpair(kk,1)) = newRassoc(tpair(kk,1),tpair(kk,2));
                end

                swappedRs = zeros(Nperm_sz,numRs);
                swappedTypes = zeros(Nperm_sz,nMax);  
                type_holder = clusters(ii,1:nMax);
                for ll = 1:Nperm_sz
                    for kk = 1:tpair_sz
                        BodyA = Nperm(ll,tpair(kk,1));
                        BodyB = Nperm(ll,tpair(kk,2));
                        swappedRs(ll,kk) = newRassoc(BodyA,BodyB);
                        if qq ~= nMax
                            swappedTypes(ll,:) = [type_holder(Nperm(ll,:)) zeros(1,nMax-qq)];
                        else
                            swappedTypes(ll,:) = type_holder(Nperm(ll,:));
                        end            
                    end
                end           
                allpossperms = [swappedTypes swappedRs];
                % Check each permutation against the list of known interactions
                % (i.e. INTERACTIONS matrix). Increase the count of that
                % interaction
                flag = 0;
                for nn = 1:Nperm_sz
                    for mm = 1:INT_sz            
                        checkdif = abs(allpossperms(nn,:) - INTERACTIONS(mm,:));
                        if all(checkdif < 0.05)             % Adjust this tolerance if using substantially large Rmax.
                            flag = 1;
                            counts(jj,mm) = counts(jj,mm) + 1;
                            break
                        elseif mm == INT_sz && nn == Nperm_sz && flag == 0
                            fprintf('A structure could not be assigned to an interaction!\n')
                        end
                    end
                    if flag == 1
                        break
                    end
                end   
            end 
        end
        fprintf('...Done!\n')
        Int_sz = size(INTERACTIONS,1);
        output_stuff = [(1:Int_sz)' INTERACTIONS];
        surf_normalized_counts(jj,:) = counts(jj,:)./Ns;
        surf_energy(jj) = energy(jj)./Ns;
        ads_energy(jj) = energy(jj)./num_adsorbates;
        % Put "analyzed" at end of the file for this configuration
        backhome = cd('./new_configs');
        fID = fopen(all_configs(jj+2).name,'a');
        fprintf(fID,'\nanalyzed\n');
        if prob_flag == 0
            fprintf(fID,'Confidence\n');
        else
            fprintf(fID,'No Confidence\n');
        end
        dlmwrite(all_configs(jj+2).name,1:size(counts,2),'delimiter','\t','-append');
        dlmwrite(all_configs(jj+2).name,counts(jj,:),'delimiter','\t','-append');
        fprintf(fID,'Coverage:     %6.5g ML\nNumber of sites: %g\nSurface Energy: %7.4g eV/site\nAdsorption Energy: %8.4',coverages(jj),Ns,surf_energy(jj),ads_energy(jj));
        fprintf(fID,'\n\nHere are the interactions types for easy reference:\n');
        dlmwrite(all_configs(jj+2).name,output_stuff,'delimiter','\t','precision','%4.3g','-append');
        fclose(fID);
        cd(backhome)
        if prob_flag == 1 && move_input3 == 1
            movefile(['./new_configs/' all_configs(jj+2).name],'./prob_configs/')
            fprintf('File Moved to "prob_configs"\n');
        elseif prob_flag == 0 && move_input2 == 1
            movefile(['./new_configs/' all_configs(jj+2).name],'./configs/')
            fprintf('File Moved to "configs"\n');
        end
    end
    if num_configs > 0
        % Extract zero coverage energy from 'zero_energy.txt'
        fID = fopen('zero_energy.txt');                           
        tline = fgetl(fID);  
        zero_en = cell2mat(textscan(tline, '%f'));
        fclose(fID);

        if no_LG == 0
            % Extract LG EICs from LG_EICs.txt
            fID = fopen('LG_EICs.txt');                           
            tline = fgetl(fID);
            line=0;
            while ischar(tline)
                line = line + 1 ;
                tline = fgetl(fID);         
            end
            line = line - 1;

            fclose(fID);
            LG_EICs = zeros(line-1,1); 
            fID = fopen('LG_EICs.txt');
            tline = fgetl(fID);
            CE = cell2mat(textscan(tline, '%f'));
            for ii = 1:line-1
                tline = fgetl(fID);
                LG_EICs(ii) = cell2mat(textscan(tline, '%f'));
            end
            fclose(fID);

            
           %%% Determine if any of the new structures are not unique
           %%% amongst themselves. Remove all not unique structures from
           %%% calculation of external CV score, keeping the lowest energy
           %%% version
            repeat_flag = 0;
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
                if totoss > 0
                    repeat_flag = 1;
                end
            end
            tokeep = setdiff(indMat,totoss);    
            surf_normalized_counts = surf_normalized_counts(tokeep,:);               
            coverages = coverages(tokeep);                                         % grab the associated coverages
            surf_energy = surf_energy(tokeep);
            newest_configs = newest_configs(tokeep);
            
            %%% Determine if any of these new structures are equivalent to
            %%% previous structures: remove these from evaluation of
            %%% external CV score
            
            if ~isempty(MAT)
               [~,repi] = intersect(surf_normalized_counts,MAT,'rows','stable'); 
               if ~isempty(repi)
                   repeat_flag = 1;
                   surf_normalized_counts(repi,:) = [];
                   relevant_counts = surf_normalized_counts(:,CE);
                   surf_energy(repi,:) = [];
                   coverages(repi) = [];
                   surf_energy = surf_energy - zero_en;
                   labels_tokeep = setdiff(1:numel(newest_configs),repi);
                   newest_configs = newest_configs(labels_tokeep);
               else
                   relevant_counts = surf_normalized_counts(:,CE);
                   surf_energy = surf_energy - zero_en;
               end
            end
            %%%%
             % Subtract off the mu_elec energy provided by the user (if it exists)
            if exist('mu_elec','var')==1
                if ~isempty(mu_elec)
                    sz_mu_elec = size(mu_elec,1);

                    for ii = 1:sz_mu_elec
                        cluster_num = mu_elec(ii,1);
                        cluster_ECI = mu_elec(ii,2);

                        part_en = surf_normalized_counts(:,cluster_num)*cluster_ECI;
                        surf_energy = surf_energy - part_en;
                    end
                end
            end
          
            if no_mf == 1
                % Calculate predicted surface energy based on this LG model
                predicted_energy = relevant_counts*LG_EICs;
            else
                % Extract mean field coefficients from mf_coeffs.txt
                ID = fopen('mf_coeffs.txt');                           
                tline = fgetl(fID);
                line=0;
                while ischar(tline)
                    line = line + 1 ;
                    tline = fgetl(fID);         
                end

                fclose(fID);
                mf_coeffs = zeros(line,1); 
                fID = fopen('mf_coeffs.txt');        
                for ii = 1:line
                    tline = fgetl(fID);
                    mf_coeffs(ii) = cell2mat(textscan(tline, '%f'));
                end
                fclose(ID);

                % Calculate the predicted mean field energy: 
                mf_energy = 0;
                for ii = 1:line
                    mf_energy = mf_energy + mf_coeffs(ii)*coverages.^ii;
                end
                % ...and then predicted surface energy:
                predicted_energy = mf_energy + relevant_counts*LG_EICs;
            end
            residuals = surf_energy - predicted_energy;
            predicted_ads_en = predicted_energy./coverages;
            ads_en = surf_energy./coverages;       
            sqrd_residuals = residuals.^2;
            ext_CV = sqrt(mean(sqrd_residuals));
            resid_stdev = sqrt(mean(sqrd_residuals)-(mean(residuals))^2);
            Rr = sqrd_residuals;
            ext_stdev = (mean(Rr.^2)-(mean(Rr))^2)^(1/4);
            if num_configs > 0
                display_flag = 1;
            end
            % Plot up the results   
            mf_model = 0;
            th = 0:0.01:1;
            zeroline = zeros(size(th,2),2);

            subplot(2,1,1)
            scatter(coverages,ads_en,'sk')
            hold on
            scatter(coverages,predicted_ads_en,'r+')
            ylabel('Ads. En. (eV/adsorbate)')
            title('Current LG Model')
            xlim([0 1])
            hold off

            subplot(2,1,2)
            scatter(coverages,residuals,'r+')
            hold on
            plot(th,zeroline,'k')
            xlabel('OH Coverage (ML)')
            ylabel('Residuals (eV/site)')
            xlim([0 1])
            hold off
        end
    end
end


% Loop over configurations in directory "configs", which should be located
% in your working directory where this script is located

all_configs = dir('./configs');         % command 'dir' assigns each file a fileID starting from 3
num_configs = numel(all_configs)-2;     % since fileID 1 and 2 are STDIN and STDERR we want to discount these
counts = zeros(num_configs,INT_sz);       % A whole bunch of preallocation follows
energy = zeros(num_configs,1);
coverages = zeros(num_configs,1);
surf_normalized_counts = counts;
surf_energy = energy;
confnum = 0;
for jj=1:num_configs                    % start working through the configurations...
    fprintf('CONFIGURATION #%d, %s:\n',jj, all_configs(jj+2).name)
    backhome = cd('./configs/');
    confnum = confnum + 1;
    analyzed_flag = 0;
    % Determine how many lines there are in this file
    fID = fopen(all_configs(jj+2).name);
    tline = fgetl(fID);
    line=0;
    emptylines = 0;
    while ischar(tline)
        if isempty(tline)
            tline = fgetl(fID);
            continue
        end

        if strcmp(tline,'analyzed')                         % check if this configuration has already been analyzed
            analyzed_flag = 1;                              % Turn the analyzed_flag on
            tline = fgetl(fID);
            if strcmp(tline,'Confidence') || strcmp(tline,'No Confidence')
                tline = fgetl(fID); %#ok<NASGU>
            end
            tline = fgetl(fID);                             % advance to where the counts should be               
            counts(jj,:) = cell2mat(textscan(tline, '%f')); % grab these counts for the matrix
            break
        else
        line = line + 1 ;
        tline = fgetl(fID);   
        end
    end   
    fclose(fID);    
     
    line = line -1;
    CELL = zeros(2,2);
    SIG = zeros(line-2,3);

    % Now extract the unit cell of the configuration (CELL) along with the
    % positions of the adsorbates (SIG) and the configuration energy
    % (energy)
    fID = fopen(all_configs(jj+2).name);
    
    for linenum = 1:2
        tline = fgetl(fID);
        CELL(linenum,:) = cell2mat(textscan(tline, '%f'));
    end
    for linenum = 3:line
        tline = fgetl(fID);
        SIG(linenum-2,:) = cell2mat(textscan(tline, '%f'));
    end
    tline = fgetl(fID);
    energy(jj) = str2double(tline);         
   
    fclose(fID);
    cd('../')


    Ns = norm(det([CELL(1,:);CELL(2,:)]))*sites_per_cell;
    rel_SIG = [];
    for ii = 1:max(size(coverage_sites))
        rel_SIG = [rel_SIG; SIG(SIG(:,3)==coverage_sites(ii),:)];  %#ok<AGROW>
    end
    num_adsorbates = size(rel_SIG,1);

    coverages(jj) = num_adsorbates/Ns;
    
    num_adsorbates = size(SIG,1);
    
    % find length of all effective vector lengths in the unit cell
    cell_length = zeros(1,4);
    cell_length(1) = CELL(1,:)*normR(1:2,1:2)*CELL(1,:)';
    cell_length(2) = CELL(2,:)*normR(1:2,1:2)*CELL(2,:)';
    cell_length(3) = (CELL(1,:)-CELL(2,:))*normR(1:2,1:2)*(CELL(1,:)-CELL(2,:))';
    cell_length(4) = (CELL(1,:)+CELL(2,:))*normR(1:2,1:2)*(CELL(1,:)+CELL(2,:))';

    % If any of those length are equal to or less than the problematic
    % length, flag this structure
    prob_flag = 0;
    if sum(cell_length-prob_length^2 < 1) > 0.999 && sum(cell_length-prob_length^2 < 1) < 1.999
        if exist('prob_configs','dir') ~= 7
            mkdir('prob_configs');
        end
        prob_flag = 1;            
    end
        
    if analyzed_flag == 1                                   % if this configuration has already been analyzed 
        fprintf('Already analyzed!\n')
        surf_normalized_counts(jj,:) = counts(jj,:)./Ns;
        surf_energy(jj) = energy(jj)./Ns;
        if prob_flag == 1
            movefile(['./configs/' all_configs(jj+2).name],'./prob_configs/')
            fprintf('File Moved to "prob_configs"\n');
        end
        continue                                            % ...then we can skip the rest of this iteration
    end
    
    % Initial Size of sig matrix (dynamic preallocation of memory)
    BLOCK = 100;
    list_size = BLOCK;
    col_BLOCK = 4;
    sig = zeros(BLOCK,col_BLOCK);
    
    % Find SigMaxX and SigMaxY: the dimensions needed to create a surface
    % large enough to encompass the maxR distance for each adsorbate
    
    diagonal_vec = (CELL(1,:)+CELL(2,:));
    SigMaxX = diagonal_vec(1)+maxX;
    SigMaxY = diagonal_vec(2)+maxY;
    
    cellvecX = CELL'\[2*SigMaxX ; -2*SigMaxY];
    cellvecY = CELL'\[-2*SigMaxX ; 2*SigMaxY];
    
    xMin = cellvecY(1);
    xMax = cellvecX(1);

    if xMax < xMin
        dumX = xMax;
        xMax = xMin;
        xMin = dumX;
    end

    xMin = floor(xMin);
    xMax = ceil(xMax);

    yMin = cellvecX(2);
    yMax = cellvecY(2);

    if yMax < yMin
        dumY = yMax;
        yMax = yMin;
        yMin = dumY;
    end

    yMin = floor(yMin);
    yMax = ceil(yMax);
    
    xTot = xMax - xMin;
    yTot = yMax - yMin;
        
    % Populate "sig" matrix
    kk = 1;
    SIG_dif = zeros(num_adsorbates,2);
    for ii = 1:num_adsorbates     
        for xx = 0:xTot
            sigx = xx + xMin;         
            for yy = 0:yTot    
                sigy = yy + yMin;
                repeat = SIG(ii,1:2)+(CELL'*[sigx;sigy])';
                SIG_dif = SIG(:,1:2) - repeat;                      % difference between this repeated position and all adsorbates' positions                   
                R_check = (diag([SIG_dif zeros(num_adsorbates,1)]*normR*[SIG_dif zeros(num_adsorbates,1)]'));               
                if any(R_check <= maxRmax+1)
                    sig_type = Site(4,SIG(ii,3));
                    Shift = Site(1:3,SIG(ii,3))';
                    sig(kk,:) = [[repeat 0]+Shift sig_type];
                    kk = kk + 1;                      
                end
                % Add new block of memory to sig matrix if needed
                check = sig(any(sig~=0,2),:);
                if size(check,1)/size(sig,1) > 0.95                 % Only 5% of the current allocation is left
                    list_size = size(check,1) + BLOCK;
                    sig(kk+1:list_size,:)=0;                        % ...so add a new block of memory
                end                         
            end
        end  
    end

    sig = sig(any(sig~=0,2),:);
    sig = sortrows(sig,1:3);

    nTot = size(sig,1);
    X = sig(:,1);
    Y = sig(:,2);
    Z = sig(:,3);
    % Find where each site within the unit cell is within the sig matrix 
    sigPos = zeros(1,num_adsorbates);
    for ii = 1:num_adsorbates
        Shift = Site(1:3,SIG(ii,3))';
        sig_type = Site(4,SIG(ii,3));
        check = [[SIG(ii,1:2) 0]+Shift sig_type];
        [~,sigPos(ii)] = intersect(sig,check,'rows');
    end       

    % All combinations of pairs for every adsorbate on the SigMaxX by
    % SigMaxY surface
    CombTot = sortrows(nchoosek(1:nTot,2),2);
    nCombs = size(CombTot,1);
    frstn = CombTot(:,1)';
    scndn = CombTot(:,2)';

    % Find the delta X,Y,Z between each pair of bodies with each site
    delX = zeros(1,nCombs);
    delY = zeros(1,nCombs);
    delZ = zeros(1,nCombs);
    delXYZ = zeros(3,nCombs);
    R = zeros(1,nCombs);
    for ii = 1:nCombs
        delX(ii) = X(scndn(ii))-X(frstn(ii));
        delY(ii) = Y(scndn(ii))-Y(frstn(ii));
        delZ(ii) = Z(scndn(ii))-Z(frstn(ii));
        delXYZ(:,ii) = [delX(ii);delY(ii);delZ(ii)];
        R2 = delXYZ(:,ii)'*normR*delXYZ(:,ii);
        R(ii) = (R2);
    end

    % Associate R's with each of its associated pair's positions.
    % Swaping these positions changes nothing so do that now.
    Rassoc = zeros(nTot);
    for h = 1:size(CombTot,1)
        Rassoc(CombTot(h,1),CombTot(h,2)) = R(h);
        Rassoc(CombTot(h,2),CombTot(h,1)) = Rassoc(CombTot(h,1),CombTot(h,2));
    end

    % Now find all pairs that contain the number of the position of each
    % supercell adsorbate. Use this to extract all the R's that connect the
    % various n-body interactions    
    clusters = zeros(BLOCK, nMax+numRs);   
    Int_sz = BLOCK;
    int_ptr = 1;

    % Create point (V naught) clusters first
    principle_sites = unique(Site(4,:));
    if size(principle_sites,1) > 1
        uu = unique(linked);
        nn = histcounts(linked);
        over_linked = uu(nn>1);
        principle_site = intersect(principle_sites,uu);

        for ii = 1:size(over_linked,2)
            over_linked_sites = linked(linked==over_linked(ii),:);
            principle_sites = [principle_sites min(over_linked_sites(:))]; %#ok<AGROW>
        end
        non_linked_sites = principle_sites;
        for ii = 1:size(linked,1)
            if isempty(intersect(linked(ii,:),principle_sites))
                non_linked_sites = [non_linked_sites min(linked(ii,:))]; %#ok<AGROW>
            end
        end
    end

        non_linked_sites = principle_sites;
    for ii = 1:size(non_linked_sites,2)
        num_sites_here = size(SIG(SIG(:,3)==ii,:),1);
        for  ww= 1:num_sites_here
            clusters(int_ptr,1) = non_linked_sites(ii);
            posclusters(int_ptr,1) = non_linked_sites(ii);
            int_ptr = int_ptr + 1;
        end
    end
    
    % Now go through and find the rest of the clusters
    for ii = 2:nMax
        fprintf('\tFinding all %d-body clusters...\n',ii)
        for kk = 1:num_adsorbates
            fprintf('\t\tSite #%d\n',kk)
            curr_ads_type = sig(sigPos(kk),4);
            vec_dif = sig(:,1:3) - sig(sigPos(kk),1:3);
            R_vec_dif = (diag(vec_dif*normR*vec_dif'));       
            sigwR = [sig R_vec_dif];
            irrel_sigA = sigwR(sig(:,1)==sig(sigPos(kk),1) & sig(:,2)<sig(sigPos(kk),2),:);       % Find the "irrelevant" positions in sig (A) those where Y < X when X = the X of adsorbate kk
            irrel_sigB = sigwR(sig(:,1)<sig(sigPos(kk),1),:);                                     % ...(B) those where X is lower than the adsorbate        
            irrel_sigD = sigwR(sig(:,3)<sig(sigPos(kk),3),:);                                     % ...(D) whose site number is less than the adsorbates
            ind_sigC = zeros(10000,1);
            ind_ptr = 1; 
            for bb = 1:sMax                                                            % ...(C) those where their distance from the adsorbate is greater than this type's n-body Rmax                           
                nd = ndims(Rmax);                                                       % Since it can change, how man dimension in the Rmax array?
                hold_ind = repmat({':'},1,nd-3);                                        % Create nd-3 ":" indices to insert into the Rmax indexing
                subsetR = Rmax(curr_ads_type,bb,hold_ind{:},ii);                        % Get the 2-body subset of Rmax of the current site type and bb
                maxSubsetR = max(subsetR(:));                
                temp_ind = find(sigwR(:,4)==bb & sigwR(:,5) > maxSubsetR); 
                sz_temp = size(temp_ind,1);
                ind_sigC(ind_ptr:(ind_ptr+sz_temp-1)) = temp_ind;
                ind_ptr = ind_ptr + sz_temp;
            end          
            ind_sigC = ind_sigC(ind_sigC ~= 0,:);
            ind_sigC = sort(ind_sigC);
            irrel_sigC = sigwR(ind_sigC,:);                                                                                                               
            irrel_sig = [irrel_sigA ; irrel_sigB ; irrel_sigC; irrel_sigD];
            irrel_sig = setdiff(irrel_sig,sigwR(sigPos(kk),:),'rows');
            [rel_sig,rel_bodies] = setdiff(sigwR,irrel_sig,'rows');                               % Find which bodies (indices in sig) fit the above criteria
            if size(rel_bodies,1) <= 1
                    continue
            end
            rel_combs = nchoosek(rel_bodies',ii);                                               % Use these bodies to find the relevant n-body combos
            rel_combs = rel_combs(any(rel_combs == sigPos(kk),2),:);                            % Restrict the relevant bodies to those with the kk adsorbate in it
            bTot = size(rel_combs,1);
            dist = zeros(bTot,numRs);
            types = reshape(sig(rel_combs',4),ii,[])';                                          % rel_bodies are read row THEN column by sig, so transpose rel_bodies first, then read their types (col 4 of sig). Reshape this back into a matrix the same shape as rel_bodies
            good_types = zeros(bTot,ii);
            zz = 1;
            for ll = 1:bTot
                pairs = sortrows(nchoosek(rel_combs(ll,:),2),2);
                num_pairs = size(pairs,1);   
                dist_temp = diag(Rassoc(pairs(:,1),pairs(:,2)))';
                index_vec = num2cell([types(ll,:) ones(1,nMax-ii) ii]);
                if all(dist_temp <= Rmax(sub2ind(size(Rmax),index_vec{:})))   
                    dist(ll,1:num_pairs) = dist_temp; 
                    good_types(ll,:) = types(ll,:);  
                end
            end                   
            fill_zeros = zeros(bTot,nMax-ii);
            dum_types = [good_types fill_zeros];              
            dum = [dum_types dist];  
            dum(all(dum == 0,2),:) = [];
            dum = dum(all(dum(:,nMax+1:nMax+num_pairs) >= Rmin-0.01,2),:);
            dum_sz = size(dum,1);            
                % Fill memory as needed
                check = clusters(all(clusters == 0,2),:);
                check_sz = size(check,1);               
                while check_sz < dum_sz
                    Int_sz = Int_sz + BLOCK;
                    clusters(int_ptr+1:Int_sz,:)=0;
                    
                    check = clusters(all(clusters == 0,2),:);
                    check_sz = size(check,1);               
                end     

            clusters(int_ptr:int_ptr+dum_sz-1,:) = dum;
            int_ptr = int_ptr + dum_sz + 1;                     
        end 

    end 
    % Clean up clusters matrix: remove unused space
    clusters = clusters(any(clusters~=0,2),:);
    TotInts = size(clusters,1);
    
    % Find where each n-body interaction starts and ends
    section_ctr = zeros(1,nMax);
    dummy_int = clusters;
    for mm = 2:nMax
        if mm ~= nMax
            section = dummy_int(all(dummy_int(section_ctr(mm-1)+1:end,mm+1)==0,2),:);
        else
            section = dummy_int(section_ctr(mm-1)+1:end,:);
        end 
        section_ctr(mm) = section_ctr(mm-1) + size(section,1);
    end
    section_ctr(nMax+1) = TotInts;
    fprintf('\nAll n-body clusters found!\n\nNow assigning them to the provided interactions...\n')    
    
    % Now permute the n-bodies and check against the INTERACTIONS matrix
    
    for qq = 2:nMax
        Nperm = sortrows(perms(1:qq));        % Find all permutations of qq bodies
        Nperm_sz = size(Nperm,1);
        tpair = sortrows(nchoosek(1:qq,2),2);              % Total combinations of pairs used here according to number qq
        tpair_sz = size(tpair,1);
        for ii = section_ctr(qq-1)+1:section_ctr(qq)
            newRassoc = zeros(nMax);
            for kk = 1:tpair_sz
                newRassoc(tpair(kk,1),tpair(kk,2)) = clusters(ii,nMax+kk);
                newRassoc(tpair(kk,2),tpair(kk,1)) = newRassoc(tpair(kk,1),tpair(kk,2));
            end

            swappedRs = zeros(Nperm_sz,numRs);
            swappedTypes = zeros(Nperm_sz,nMax);  
            type_holder = clusters(ii,1:nMax);
            for ll = 1:Nperm_sz
                for kk = 1:tpair_sz
                    BodyA = Nperm(ll,tpair(kk,1));
                    BodyB = Nperm(ll,tpair(kk,2));
                    swappedRs(ll,kk) = newRassoc(BodyA,BodyB);
                    if qq ~= nMax
                        swappedTypes(ll,:) = [type_holder(Nperm(ll,:)) zeros(1,nMax-qq)];
                    else
                        swappedTypes(ll,:) = type_holder(Nperm(ll,:));
                    end            
                end
            end           
            allpossperms = [swappedTypes swappedRs];
            % Check each permutation against the list of known interactions
            % (i.e. INTERACTIONS matrix). Increase the count of that
            % interaction
            flag = 0;
            for nn = 1:Nperm_sz
                for mm = 1:INT_sz            
                    checkdif = abs(allpossperms(nn,:) - INTERACTIONS(mm,:));
                    if all(checkdif < 0.05)             % Adjust this tolerance if using substantially large Rmax.
                        flag = 1;
                        counts(jj,mm) = counts(jj,mm) + 1;
                        break
                    elseif mm == INT_sz && nn == Nperm_sz && flag == 0
                        fprintf('A structure could not be assigned to an interaction!\n')
                    end
                end
                if flag == 1
                    break
                end
            end   
        end 
    end
    fprintf('...Done!\n')
    Int_sz = size(INTERACTIONS,1);
    output_stuff = [(1:Int_sz)' INTERACTIONS];
    surf_normalized_counts(jj,:) = counts(jj,:)./Ns;
    surf_energy(jj) = energy(jj)./Ns;
    ads_energy(jj) = energy(jj)./num_adsorbates;
        % Put "analyzed" at end of the file for this configuration
        backhome = cd('./configs');
        fID = fopen(all_configs(jj+2).name,'a');
        fprintf(fID,'\nanalyzed\n');
        if prob_flag == 0
            fprintf(fID,'Confidence\n');
        else
            fprintf(fID,'No Confidence\n');
        end
        dlmwrite(all_configs(jj+2).name,1:size(counts,2),'delimiter','\t','-append');
        dlmwrite(all_configs(jj+2).name,counts(jj,:),'delimiter','\t','-append');
        fprintf(fID,'Coverage:     %6.5g ML\nNumber of sites: %g\nSurface Energy: %7.4g eV/site\nAdsorption Energy: %8.4',coverages(jj),Ns,surf_energy(jj),ads_energy(jj));
    fprintf(fID,'\n\nHere are the interactions types for easy reference:\n');
    dlmwrite(all_configs(jj+2).name,output_stuff,'delimiter','\t','precision','%4.3g','-append');
    fclose(fID);
    cd(backhome)
end
% Subtract off the zero coverage (i.e. clean surface) energy 
% If this energy is not available warn the user
no_zero_flag = 0;
if any(ismember(coverages,0,'rows')) == 1
    [~,iz] = intersect(coverages,0);
    zero_en = surf_energy(iz);
    surf_energy = surf_energy - zero_en;
    
    fID = fopen('zero_energy.txt','w');
    fprintf(fID,'%12.8f',zero_en);
    fclose(fID);
else
    no_zero_flag = 1;
end

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
surf_normalized_counts = surf_normalized_counts(tokeep,:);
coverages = coverages(tokeep);                                         % grab the associated coverages
surf_energy = surf_energy(tokeep);


atat_names = all_configs(3:end);
atat_names = {atat_names.name};
atat_names = atat_names(tokeep);
% Subtract off the mu_elec energy provided by the user (if it exists)
part_en = 0;
if exist('mu_elec','var')==1
    if ~isempty(mu_elec)
        sz_mu_elec = size(mu_elec,1);

        for ii = 1:sz_mu_elec
            cluster_num = mu_elec(ii,1);
            cluster_ECI = mu_elec(ii,2);

            part_en = part_en + surf_normalized_counts(:,cluster_num)*cluster_ECI;
            surf_energy = surf_energy - part_en;
        end
    end
end
ads_normalized_counts = surf_normalized_counts./coverages;
ads_normalized_counts(~isfinite(ads_normalized_counts(:))) = 0;
ads_energy = surf_energy./coverages;
ads_energy(~isfinite(ads_energy)) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% v7 addition %%%%%%%%%%%%%%%%%%%

% Find isosteric ground states
[uniq_covs] = unique(round(coverages*1E9)/1E9);             % i.e. unique out to 9 decimal points
num_uniq_covs = numel(uniq_covs);
iso_gs = (1:num_uniq_covs)*0;


for ii = 1:num_uniq_covs
    subset = find(abs(coverages - uniq_covs(ii)) <= 1E-8);     % grab indices of coverages equal to the ii'th unique coverage        
    subset_energies = surf_energy(subset);                     % Find the energies for their corresponding structures
    [~,isu] = min(subset_energies);                                 % which index within subset energies corresponds to the lowest energy structure (the isosteric ground state)
    iso_gs(ii) = subset(isu);                                  % and this corresponds to which index within the entire coverages vector
end

iso_gs_ens = surf_energy(iso_gs);                                   % grab the energies of the iso_gs
iso_gs_covs = coverages(iso_gs);                                    % grab the coverages of the iso_gs


% Find true ground states
is_true_gs = (1:num_uniq_covs)*0;                                   % preallocating (is_true_gs = "is this a true ground state")
is_true_gs(1) = 1;                                                  % the first (empty coverage) iso_gs is always a true_gs

ii = 1;
while ii < num_uniq_covs
    Slope = (1:num_uniq_covs)*inf;                                  % preallocating
    for jj = ii+1:num_uniq_covs                                     % look at all slopes connecting this iso_gs to all future iso_gs
        Slope(jj) = (iso_gs_ens(jj) - iso_gs_ens(ii))/(iso_gs_covs(jj) - iso_gs_covs(ii));    % slope (a forward finite difference)
    end
    minS = min(Slope);                                              % the minimum slope
    minjj = find(Slope == minS,1);                                        % find which iso_gs is the next true_gs
    is_true_gs(minjj) = 1;                                          % mark this is as a true_gs
    
    ii = minjj;                                                     % move 'ii' counter to the location of the new true_gs
end

true_gs = iso_gs(is_true_gs == 1);                                  % find the indices of the true ground states

%%% Now that we know the final surf energies we can grab those
%%% corresponding to the ground states %%%

true_gs_covs = coverages(true_gs);                                  % find the corresponding true_gs coverages
true_gs_ens = surf_energy(true_gs);                                 % find the corresponding true_gs surf energies
true_gs_ads_en = ads_energy(true_gs);                               % find the corresponding true_gs adsorption energies
true_gs_names = atat_names(true_gs);                                % grab the true ground states' structure names

true_gs_slopes = true_gs.*0;                                        % preallocating for the gs_slopes
for ii = 2:numel(true_gs)
    true_gs_slopes(ii) = (true_gs_ens(ii) - true_gs_ens(ii-1))/(true_gs_covs(ii) - true_gs_covs(ii-1)); %forward finite difference slope
end

form_E = surf_energy - surf_energy(coverages == max(coverages))*coverages/max(coverages);   % calculate the formaiton energy

convex_hull = true_gs_ens - surf_energy(coverages == max(coverages))*true_gs_covs/max(coverages);   % do the same for the ground states to get convex hull

convex_hull_ads_E = true_gs_ens./true_gs_covs;                      % adsorption energy convex hull
%convex_hull_ads_E(~isfinite(convex_hull_ads_E)) = [];               % convert NaN due to division by zero to "0"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v8 Find the normal matrix of constraints
if ads_flag == 1
    temp_counts = ads_normalized_counts;
    temp_energy = ads_energy;
else
    temp_counts = surf_normalized_counts;
    temp_energy = surf_energy;
end

% First all structures except the ground states
constraints = temp_counts;                       % preallocate
delta_else = temp_energy;                            % preallocate
kk = 2;
for ii = 1:num_uniq_covs 
    phi = (uniq_covs(ii) - true_gs_covs(kk-1))/(true_gs_covs(kk) - true_gs_covs(kk-1));         % fractional distance between the previous gs and the next
    if phi - 1 > 0                                                  % if phi is greater than 1, we are between two new ground states 
        kk = kk + 1;                                                % ...and must increment to that pair
        phi = (uniq_covs(ii) - true_gs_covs(kk-1))/(true_gs_covs(kk) - true_gs_covs(kk-1));         % fractional distance between the previous gs and the next
    end
   
    subset = find(abs(coverages - uniq_covs(ii)) <= 1E-8);      % grab indices of coverages equal to the ii'th unique coverage        
    subset_energies = temp_energy(subset);                      % Find the energies for their corresponding structures
    subset_counts = temp_counts(subset,:);             % the counts for this isosteric subset
    
    if abs(phi - 1) < 1E-8
        gs_line = temp_counts(kk,:);               % if phi = 1 then we're at ground state kk and can just use the ground state directly
        delta_else(subset) = 0;
    else
        gs_line = temp_counts(kk-1,:) + phi*(temp_counts(kk,:) - temp_counts(kk-1,:));  % counts that produce the predicted point on the ground state line
        delta_else(subset) = temp_energy(iso_gs(ii)) - (temp_energy(true_gs(kk-1)) + phi*(temp_energy(true_gs(kk)) - temp_energy(true_gs(kk-1)))); % use the actual distance between iso_gs and the gs line as the minimum constrained distance
    end
    constraints(subset,:) = subset_counts - gs_line;
end
% And then the ground states
constraints(true_gs,:) = temp_counts(true_gs,:);
delta_else(true_gs) = 0;
energy_gs = temp_energy(true_gs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now print out a summary and plot
fprintf('----------------------------------------------------------\n')
fprintf('The provided dataset exhibits the following ground states:\n')
T = table(true_gs_names',true_gs_covs,true_gs_ads_en, true_gs_ens,true_gs_slopes');
T.Properties.VariableNames = {'Structure' 'Coverage' 'Ads_Energy' 'Surf_Energy','GS_Chem_Pot'};
disp(T)

figure
plot(coverages(coverages~=0), ads_energy(coverages~=0), 'rx')
hold on
plot(true_gs_covs, convex_hull_ads_E, 'b-o')
hold off
title('Adsorption Energy')

figure
plot(coverages, form_E, 'rx')
hold on
plot(true_gs_covs, convex_hull, 'b-o')
hold off
title('Formation Energy')

figure
plot(coverages,surf_energy,'rx')
hold on
plot(true_gs_covs,true_gs_ens,'b-o')
hold off
title('Surface Energy')



% Subtract off energy due to specified ECI values (if they exist)

part_en = 0;
part_ads_en = 0;
part_gs_slope_en = 0;
part_slope_en = 0;
if exist('ECI_val','var')==1
    if ~isempty(ECI_val)
        sz_ECI_vals = size(ECI_val,1);

        for ii = 1:sz_ECI_vals
            cluster_num = ECI_val(ii,1);
            cluster_ECI = ECI_val(ii,2);

            part_en = part_en + surf_normalized_counts(:,cluster_num)*cluster_ECI;
            part_ads_en = part_ads_en + ads_normalized_counts(:,cluster_num)*cluster_ECI;
            
                orig_surf_normalized_counts = surf_normalized_counts;   % Don't think I need this anymore...
            surf_normalized_counts(:,cluster_num) = 0; 
            ads_normalized_counts(:,cluster_num) = 0;
            
        end
    end
end
surf_energy = surf_energy - part_en;
ads_energy = ads_energy - part_ads_en;
ads_energy(~isfinite(ads_energy)) = 0;

dlmwrite('COUNTS_MATRIX.txt',counts,'delimiter',' ');
dlmwrite('NORMALIZED_MATRIX.txt',surf_normalized_counts,'delimiter',' ');
dlmwrite('COVERAGES.txt',coverages,'delimiter',' ');
dlmwrite('SURF_ENERGY',surf_energy,'delimiter',' ');

fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
fprintf('FINSIHED!\nEach configuration file has had its results added to it.\nThe strict total counts have been written to "COUNTS_MATRIX.txt"\nThe normalized counts have been written to "NORMALIZED_MATRIX.txt"\nThe coverage of each system has been written to "COVERAGES.txt"\nThe surface normalized energies have been written to "SURF_ENERGY.txt"\n')
if no_zero_flag == 1
    fprintf('\nWARNING. You have not provided a zero-coverage structure in the configs directory.\nThe current results cannot be fit to a lattice gas model!\n\n')
end
fprintf('----------------------------------------------------------------\n')
fprintf('For ease of manual manipulation of the data,\nall the above data can be found in the following variables:\n"counts"\n"surf_normalized_counts"\n"coverages"\n"surf_energy"\n')
fprintf('----------------------------------------------------------------\n')
if exist('repeat_flag','var') == 1
    if repeat_flag == 1
        if numel(repi)+numel(totoss) == 1
            fprintf('\t\tOOPS!\n\t\tBased on the available data,\n\t\t%d of these new structures is equivalent to one of\n\t\tthe structures that have already been added!\n\t\tThis structure will be ignored in evaluation of the external CV-score\n',numel(repi)+numel(totoss))
        elseif numel(repi)+numel(totoss) > 1
            fprintf('\t\tOOPS!\n\t\tBased on the available data,\n\t\t%d of these new structures are equivalent to each other\n\t\tor to structures that have already been added!\n\t\tThese structures will be ignored in evaluation of the external CV-score\n',numel(repi)+numel(totoss))            
        end
    else
        fprintf('\t\tCongrats!\n\t\tBased on the available data,\n\t\tALL of these new structures are unique!\n')
    end
end
fprintf('----------------------------------------------------------------')
if display_flag == 1
    fprintf('\nThe external-CV score is %8.6g eV/site\n',ext_CV);
    fprintf('The standard deviation of the residuals is %8.6g eV/site\n',resid_stdev);
    fprintf('The external-CV deviation is %8.6g eV/site\n',ext_stdev);
    fprintf('The residuals are:\n');
    for ii = 1:size(residuals,1)
        fprintf('%s   %8.6g \n', newest_configs{ii}, residuals(ii));
    end
    
end
%Timing stuff
elapsed=toc;
inmin = elapsed/60;
fprintf('\nThis run took %9.2f seconds (or %5.4f minutes) to run.\n',elapsed,inmin)
fprintf('----------------------------------------------------------------')
fprintf('\n----------------------------------------------------------------\n')
clearvars -except constraints energy_gs delta_else slope_counts slope_energy part_slope_en true_gs_slopes_counts true_gs_slopes part_ads_en ads_normalized_counts true_gs_ens true_gs_ads_ens true_gs_covs true_gs iso_gs diff_ads_counts diff_surf_counts diff_surf_energies diff_ads_energies true_gs_covs true_gs_ads_e  true_gs_ens true_gs_slopes ads_energy part_en orig_surf_normalized_counts atat_names newest_configs residuals mf_resids counts surf_normalized_counts coverages surf_energy sMax INTERACTIONS

