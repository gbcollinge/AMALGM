%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% INTERACTIONS_GEN.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build and find interaction terms in a cluster expansion for adsorbates
% on metal surfaces.

% v3: Hopefully output MC_POSITIONS.txt for use in MC simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
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
% you should provide the c vector as it appears in your POSCAR divided by 
% the 1st NN distance. Otherwise, you can (and should) leave it as 
% [0 ; 0 ; 1]. Mind, whatever length a unit vector within this coordinate 
% system is will be the "natural unit" used from here on out.
% NOTE: This is the only place where cartesian coordinates should be 
% encountered!

ux = [1 ; 0 ; 0];
uy = [-1/2 ; sqrt(3)/2 ; 0];
uz = [0 ; 0 ; 1];

% Alternatively, change infile to "1" and provide a file called
% "NATURAL_COORDINATES.txt" with each ux uy and uz provided as column
% vectors. 

infile = 0;

% This file must be written in decimal (floating point) format.
% An example for FCC(11) or HCP(0001) follows:
%
%     1 0 0
%     -0.5 0.866025403784439 0
%     0 0 1    
%
%     end of example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% User specified overall maximum N-body clusters to include (even if the
% max is different for different site types, still specify the max of all
% types here)

maxNbody= 5;

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

Rmin = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%& END USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% manipulate/fill in Rmax array
for ii=2:maxNbody
    Ncombs = nmultichoosek(unique(Site(4,:)),ii);
    if size(Ncombs,1) == 1
        Ncombs = ones(1,ii)*unique(Site(4,:));
    end
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

%if size(Rmax,2) ~= max(max(Nbody)) - 1
    %error('Error. Number of defined maximum distances per n-body interaction inconsistent with the number of n-body interactions specified.');
%end
fprintf('----------------------------------------------------------------\n')
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
maxRmax = max(Rmax(:));
unitX = [1;0;0];
unitY = [0;1;0];
lengthX = unitX'*normR*unitX;
lengthY = unitY'*normR*unitY;
maxX = ceil(sqrt(maxRmax)/lengthX);
maxY = ceil(sqrt(maxRmax)/lengthY);
% Total number of sites and bodies
Site = Site(:,Site(4,:)~=0);        % get rid of Sites where the "type" is 0
sMax = size(unique(Site(4,:)),2);
nMax = maxNbody;

% Total number of body-to-body pair distances
numRs = nchoosek(nMax,2);
site_max = size(Site,2);
% The "CELL" and "SIG" matrix:
CELL = [1 0 ; 0 1 ];
SIG = [zeros(site_max,2) (1:site_max)'];

num_adsorbates = size(SIG,1);

% Initial Size of sig matrix (dynamic preallocation of memory)
BLOCK = 100;
col_BLOCK = 4;
sig = zeros(BLOCK,col_BLOCK);

% Find SigMaxX and SigMaxY: the dimensions needed to create a surface
% large enough to encompass the maxR distance for each adsorbate
SigMaxX = maxX;
SigMaxY = maxY;

cellvecX = CELL'\[2*SigMaxX ; -2*SigMaxY];
cellvecY = CELL'\[-2*SigMaxX ; 2*SigMaxY];

xMin = cellvecY(1);
xMax = cellvecX(1);

xMin = sign(xMin)*ceil(abs(xMin));
xMax = sign(xMax)*ceil(abs(xMax));

yMin = cellvecX(2);
yMax = cellvecY(2);

yMin = sign(yMin)*ceil(abs(yMin));
yMax = sign(yMax)*ceil(abs(yMax));

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
            if any(R_check <= maxRmax+0.1)
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
sig_coord = zeros(num_adsorbates,4);
for ii = 1:num_adsorbates
    Shift = Site(1:3,SIG(ii,3))';
    sig_type = Site(4,SIG(ii,3));
    check = [[SIG(ii,1:2) 0]+Shift sig_type];
    [sig_coord(ii,:),sigPos(ii)] = intersect(sig,check,'rows');
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
posclusters = zeros(BLOCK, 4*nMax+numRs);     % create a matrix with enough position for each n-body, each R, and each one's 3 coordinates
min_clusters = zeros(1,BLOCK);
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
    clusters(int_ptr,1) = non_linked_sites(ii);
    posclusters(int_ptr,1) = non_linked_sites(ii);
    int_ptr = int_ptr + 1;
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
        ind_sigC = zeros(10000,1);
        ind_ptr = 1; 
        for bb = 1:sMax                                                             % ...(C) those where their distance from the adsorbate is greater than this type's n-body Rmax                           
            nd = ndims(Rmax);                                                       % Since it can change, how many dimension in the Rmax array?
            hold_ind = repmat({':'},1,nd-3);                                        % Create nd-3 ":" indices to insert into the Rmax indexing
            subsetR = Rmax(curr_ads_type,bb,hold_ind{:},ii);                        % Get the 2-body subset of Rmax of the current site type and bb
            maxSubsetR = max(subsetR(:));                                           % Find the maximum value for this subset
            temp_ind = find(sigwR(:,4)==bb & sigwR(:,5) > maxSubsetR);              % Grab the sigwR positions that are outside the maxRmax distance from the current position
            sz_temp = size(temp_ind,1);                                             % Number of these positions?
            ind_sigC(ind_ptr:(ind_ptr+sz_temp-1)) = temp_ind;                       % Add their indices to ind_sigC
            ind_ptr = ind_ptr + sz_temp;                                            % advance pointer
        end          
        ind_sigC = ind_sigC(ind_sigC ~= 0,:);
        ind_sigC = sort(ind_sigC);
        irrel_sigC = sigwR(ind_sigC,:);                                                                                                               
        irrel_sig = [irrel_sigA ; irrel_sigB ; irrel_sigC];
        irrel_sig = setdiff(irrel_sig,sigwR(sigPos(kk),:),'rows');
        [rel_sig,rel_bodies] = setdiff(sigwR,irrel_sig,'rows');                               % Find which bodies (indices in sig) fit the above criteria
        rel_combs = nchoosek(rel_bodies',ii);                                               % Use these bodies to find the relevant n-body combos
        rel_combs = rel_combs(any(rel_combs == sigPos(kk),2),:);                            % Restrict the relevant bodies to those with the kk adsorbate in it
        bTot = size(rel_combs,1);
        dist = zeros(bTot,numRs);
        types = reshape(sig(rel_combs',4),ii,[])';                                          % rel_bodies are read row THEN column by sig, so transpose rel_bodies first, then read their types (col 4 of sig). Reshape this back into a matrix the same shape as rel_bodies
        rel_pos = [];
        for yyy = 1:ii
            rel_pos = [rel_pos sig(rel_combs(:,yyy)',1:3)]; %#ok<AGROW>
        end
        good_types = zeros(bTot,ii);
        good_pos = zeros(bTot,3*ii);
        zz = 1;
        for ll = 1:bTot
            pairs = sortrows(nchoosek(rel_combs(ll,:),2),2);
            num_pairs = size(pairs,1);   
            dist_temp = diag(Rassoc(pairs(:,1),pairs(:,2)))';            
            index_vec = num2cell([types(ll,:) ones(1,nMax-ii) ii]);
            if all(dist_temp <= Rmax(sub2ind(size(Rmax),index_vec{:})))  
                dist(ll,1:num_pairs) = dist_temp; 
                good_types(ll,:) = types(ll,:);  
                good_pos(ll,:) = rel_pos(ll,:);
            end
        end                   
        fill_zeros = zeros(bTot,nMax-ii);
        fill_zeros2 = zeros(bTot, 3*(nMax-ii));
        dum_types = [good_types fill_zeros];  
        dum_pos = [good_pos fill_zeros2];
        dum = [dum_types dist];  
        dum2 = [dum dum_pos];
        dum(all(dum == 0,2),:) = [];
        dum2(all(dum2 == 0,2),:) = [];
        dum = dum(all(dum(:,nMax+1:nMax+num_pairs) >= Rmin-0.01,2),:);
        dum_sz = size(dum,1);           
            % Fill memory as needed
            check = clusters(all(clusters == 0,2),:);           
            check_sz = size(check,1);               
            while check_sz < dum_sz
                Int_sz = Int_sz + BLOCK;
                clusters(int_ptr+1:Int_sz,:)=0;
                posclusters(int_ptr+1:Int_sz,:)=0;
                min_clusters(int_ptr+1:Int_sz) =0;
                check = clusters(all(clusters == 0,2),:);
                check_sz = size(check,1);               
            end     
        
        clusters(int_ptr:int_ptr+dum_sz-1,:) = dum;
        posclusters(int_ptr:int_ptr+dum_sz-1,:) = dum2;
        int_ptr = int_ptr + dum_sz + 1;                     
    end 
    
end 
% Clean up clusters matrix: remove unused space
clusters = clusters(any(clusters~=0,2),:);
posclusters = posclusters(any(posclusters~=0,2),:);
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
fprintf('\nAll n-body clusters found!\n\nNow reducing to non-equivalent interactions...\n\n')    

% Initial Size of INTERACTIONS matrix (dynamic preallocation of memory)
BLOCK = 100;
list_size = BLOCK;
col_BLOCK = nMax+numRs;
INTERACTIONS = zeros(BLOCK,col_BLOCK);
POS_INTERACTIONS = zeros(BLOCK,4*nMax+numRs);
list_ptr = 1;

% Now permute the n-bodies and check against all previously found
% interactions, if it's a new one, add it to INTERACTIONS

for qq = 2:nMax
    fprintf('\n\nFinding equivalent site-permutations of %d-body clusters\nCounting unique clusters as they are found:\n',qq)   
    Nperm = sortrows(perms(1:qq));        % Find all permutations of qq bodies
    Nperm_sz = size(Nperm,1);
    tpair = sortrows(nchoosek(1:qq,2),2); % Total combinations of pairs used here according to number qq
    tpair_sz = size(tpair,1);
    for ii = section_ctr(qq-1)+1:section_ctr(qq)
        if size(INTERACTIONS(any(INTERACTIONS ~= 0,2),:),1) == 0
            fprintf('1...')
            INTERACTIONS(list_ptr+1,:) = clusters(1,:);
            POS_INTERACTIONS(list_ptr+1,:) = posclusters(1,:);
            list_ptr = list_ptr + 1;
            continue
        end
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
        % Check each permutation against the interactions already found
        flag = 0;
            
        check = INTERACTIONS(any(INTERACTIONS~=0,2),:);
        INT_sz = size(check,1);    
        for nn = 1:Nperm_sz
            for mm = 1:INT_sz                          
                checkdif = abs(allpossperms(nn,:) - check(mm,:));
                if all(checkdif < 0.05)
                    flag = 1;
                    break
                elseif nn == Nperm_sz && mm == INT_sz && flag ==0
                    fprintf('%d...',list_ptr)
                    if mod(list_ptr,7)==0
                        fprintf('\n')
                    end
                    INTERACTIONS(list_ptr+1,:) = clusters(ii,:);
                    POS_INTERACTIONS(list_ptr+1,:) = posclusters(ii,:);
                    list_ptr = list_ptr + 1;
                    % Add memory as it's needed
                    if list_ptr+BLOCK/20 > list_size
                        list_size = list_size + BLOCK;
                        INTERACTIONS(list_ptr+1:list_size,:)=0;
                        POS_INTERACTIONS(list_ptr+1:list_size,:)=0;
                    end
                end
            end
            if flag == 1
                break
            end
        end   
    end 
end
blah = INTERACTIONS(any(INTERACTIONS~=0,2),:);
pblah = POS_INTERACTIONS(any(POS_INTERACTIONS~=0,2),:);
blah = [blah sum(blah(:,nMax+1:end),2)];

% Create headings
Combs = sortrows(nchoosek(1:nMax,2),2);
numCombs = size(Combs,1);
for i = 1:nMax
    bawd{i} = sprintf('B%d',i); %#ok<SAGROW>
end
for i = 1:numCombs
    dists{i} = sprintf('R%d%d',Combs(i,1),Combs(i,2)); %#ok<SAGROW>
end
for i = 1:3:3*nMax
    poses{i} = sprintf('X%d',ceil(i/3)); %#ok<SAGROW>
    poses{i+1} = sprintf('Y%d',ceil(i/3)); %#ok<SAGROW>
    poses{i+2} = sprintf('Z%d',ceil(i/3)); %#ok<SAGROW>
end
    
Headings = [bawd dists];
Headings2 = [bawd dists poses];
% Get rid of unused rows
clusters(list_ptr:end,:)=[];
posclusters(list_ptr:end,:)=[];
new_inter = zeros(size(clusters,1),size(clusters,2));
pos_new_inter = zeros(size(posclusters));
tblah = blah;
ptblah = pblah;
front_pntr = 1;

% Sort the output so that interactions are clearly separated into the 
% various n-body interactions and then increase in "size" 
for i = 2:nMax
    zrow = nMax + nchoosek(i,2) + 1;
    if i == nMax
        [tempA, ti] = sortrows(tblah,[nMax:-1:1 size(tblah,2)]);
        tempB = ptblah(ti,:);
    else
        [tempA, ti] = sortrows(tblah(tblah(:,zrow)==0,:),[nMax:-1:1 size(tblah,2)]);
        tempB = ptblah(ptblah(:,zrow)==0,:);
        tempB = tempB(ti,:);
    end
    back_pntr = front_pntr+size(tempA,1)-1;
    tblah = setdiff(tblah,tempA,'rows');
    ptblah = setdiff(ptblah,tempB,'rows');
    new_inter(front_pntr:back_pntr,:) = tempA(:,1:end-1);
    pos_new_inter(front_pntr:back_pntr,:) = tempB;
    front_pntr = back_pntr + 1;
end
num_interactions = size(blah,1);

%%%% Generate MC_POSITIONS.txt %%%%
fprintf('\nWorking on MC_POSITIONS.txt now...\n\n')

clusters = zeros(BLOCK, nMax+numRs);
posclusters = zeros(BLOCK, 4*nMax+numRs+4);     % create a matrix with enough position for each n-body, each R, and each one's 3 coordinates
min_clusters = zeros(1,BLOCK);
Int_sz = BLOCK;
int_ptr = 1;

for ii = 2:nMax
    fprintf('\tFinding all %d-body clusters...\n',ii)
    for kk = 1:num_adsorbates
        fprintf('\t\tSite #%d\n',kk)
        curr_ads_type = sig(sigPos(kk),4);
        vec_dif = sig(:,1:3) - sig(sigPos(kk),1:3);
        R_vec_dif = (diag(vec_dif*normR*vec_dif'));       
        sigwR = [sig R_vec_dif];
        irrel_sigA = [];
        irrel_sigB = [];      
        ind_sigC = zeros(10000,1);
        ind_ptr = 1; 
        for bb = 1:sMax                                                             % ...(C) those where their distance from the adsorbate is greater than this type's n-body Rmax                           
            nd = ndims(Rmax);                                                       % Since it can change, how many dimension in the Rmax array?
            hold_ind = repmat({':'},1,nd-3);                                        % Create nd-3 ":" indices to insert into the Rmax indexing
            subsetR = Rmax(curr_ads_type,bb,hold_ind{:},ii);                        % Get the 2-body subset of Rmax of the current site type and bb
            maxSubsetR = max(subsetR(:));                                           % Find the maximum value for this subset
            temp_ind = find(sigwR(:,4)==bb & sigwR(:,5) > maxSubsetR);              % Grab the sigwR positions that are outside the maxRmax distance from the current position
            sz_temp = size(temp_ind,1);                                             % Number of these positions?
            ind_sigC(ind_ptr:(ind_ptr+sz_temp-1)) = temp_ind;                       % Add their indices to ind_sigC
            ind_ptr = ind_ptr + sz_temp;                                            % advance pointer
        end          
        ind_sigC = ind_sigC(ind_sigC ~= 0,:);
        ind_sigC = sort(ind_sigC);
        irrel_sigC = sigwR(ind_sigC,:);                                                                                                               
        irrel_sig = [irrel_sigA ; irrel_sigB ; irrel_sigC];
        irrel_sig = setdiff(irrel_sig,sigwR(sigPos(kk),:),'rows');
        [rel_sig,rel_bodies] = setdiff(sigwR,irrel_sig,'rows');                               % Find which bodies (indices in sig) fit the above criteria
        rel_combs = nchoosek(rel_bodies',ii);                                               % Use these bodies to find the relevant n-body combos
        rel_combs = rel_combs(any(rel_combs == sigPos(kk),2),:);                            % Restrict the relevant bodies to those with the kk adsorbate in it
        bTot = size(rel_combs,1);
        dist = zeros(bTot,numRs);
        types = reshape(sig(rel_combs',4),ii,[])';                                          % rel_bodies are read row THEN column by sig, so transpose rel_bodies first, then read their types (col 4 of sig). Reshape this back into a matrix the same shape as rel_bodies
        rel_pos = [];
        for yyy = 1:ii
            rel_pos = [rel_pos sig(rel_combs(:,yyy)',1:3)]; %#ok<AGROW>
        end
        good_types = zeros(bTot,ii);
        good_pos = zeros(bTot,3*ii);
        zz = 1;
        for ll = 1:bTot
            pairs = sortrows(nchoosek(rel_combs(ll,:),2),2);
            num_pairs = size(pairs,1);   
            dist_temp = diag(Rassoc(pairs(:,1),pairs(:,2)))';            
            index_vec = num2cell([types(ll,:) ones(1,nMax-ii) ii]);
            if all(dist_temp <= Rmax(sub2ind(size(Rmax),index_vec{:})))  
                dist(ll,1:num_pairs) = dist_temp; 
                good_types(ll,:) = types(ll,:);  
                good_pos(ll,:) = rel_pos(ll,:);
            end
        end                   
        fill_zeros = zeros(bTot,nMax-ii);
        fill_zeros2 = zeros(bTot, 3*(nMax-ii));
        dum_types = [good_types fill_zeros];  
        dum_pos = [good_pos fill_zeros2];
        dum = [dum_types dist];  
        all_sig_coord = ones(bTot,4).*sig_coord(kk,:);
        dum2 = [dum dum_pos all_sig_coord];
        dum(all(dum == 0,2),:) = [];
        dum2(all(dum2(:,1:end-4) == 0,2),:) = [];
        dum = dum(all(dum(:,nMax+1:nMax+num_pairs) >= Rmin-0.01,2),:);
        dum2 = dum2(all(dum2(:,nMax+1:nMax+num_pairs) >= Rmin-0.01,2),:);
        dum_sz = size(dum,1);
        dum_sz2 = size(dum2,1);
            % Fill memory as needed
            check = clusters(all(clusters == 0,2),:);           
            check_sz = size(check,1);               
            while check_sz < dum_sz
                Int_sz = Int_sz + BLOCK;
                clusters(int_ptr+1:Int_sz,:)=0;
                posclusters(int_ptr+1:Int_sz,:)=0;
                min_clusters(int_ptr+1:Int_sz) =0;
                check = clusters(all(clusters == 0,2),:);
                check_sz = size(check,1);               
            end     
        
        clusters(int_ptr:int_ptr+dum_sz-1,:) = dum;
        posclusters(int_ptr:int_ptr+dum_sz-1,:) = dum2;        
        int_ptr = int_ptr + dum_sz + 1;                     
    end 
    
end 
% Clean up clusters matrix: remove unused space
clusters = clusters(any(clusters~=0,2),:);
posclusters = posclusters(any(posclusters~=0,2),:);
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
fprintf('\nAll n-body clusters found!\n\n')    

% Initial Size of INTERACTIONS matrix (dynamic preallocation of memory)
BLOCK = 100;
list_size = BLOCK;
col_BLOCK = nMax+numRs;
INTERACTIONS2 = zeros(BLOCK,col_BLOCK);
POS_INTERACTIONS2 = zeros(BLOCK,4*nMax+numRs+4);
list_ptr = 1;

% Now permute the n-bodies and check against all previously found
% interactions, if it's a new one, add it to INTERACTIONS

pnt = ones(INT_sz,1);
MC_POSITIONS = zeros(100,size(posclusters,2),INT_sz); % Initialize MC_POSITIONS 3D array (1st_ind: position of NNs ; 2nd_ind: X Y Z of each NN; 1st_ind: cluster #)   

for qq = 2:nMax
    fprintf('\n\nFinding equivalent site-permutations of %d-body clusters\nCounting unique clusters as they are found:\n',qq)   
    Nperm = sortrows(perms(1:qq));        % Find all permutations of qq bodies
    Nperm_sz = size(Nperm,1);
    tpair = sortrows(nchoosek(1:qq,2),2); % Total combinations of pairs used here according to number qq
    tpair_sz = size(tpair,1);
    for ii = section_ctr(qq-1)+1:section_ctr(qq)
        if size(INTERACTIONS2(any(INTERACTIONS2 ~= 0,2),:),1) == 0

            INTERACTIONS2(list_ptr+1,:) = clusters(1,:);
            POS_INTERACTIONS2(list_ptr+1,:) = posclusters(1,:);
            list_ptr = list_ptr + 1;
            continue
        end
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
        INT_sz = size(new_inter,1);
     
        for nn = 1:Nperm_sz
            for mm = 1:INT_sz            
                checkdif = abs(allpossperms(nn,:) - new_inter(mm,:));
                if all(checkdif < 0.05)             % Adjust this tolerance if using substantially large Rmax.
                    flag = 1;
                    MC_POSITIONS(pnt(mm),:,mm) = posclusters(ii,:);
                    pnt(mm) = pnt(mm) + 1;
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
mx_pnt = max(pnt);
MC_POSITIONS(mx_pnt:end,:,:)=[];
% Fix output for use in MC simulations where coordinates are of the type 
% [X Y Z S] where S is the site type and X Y Z correspond to the INTEGER
% cell coordinate (unshifted)

sz1 = size(MC_POSITIONS,1);
sz2 = size(MC_POSITIONS,2);
sz3 = size(MC_POSITIONS,3);
bg = nMax+numRs;

for ii = 1:sz3
    for jj = 1:sz1
        dm = MC_POSITIONS(jj,:,ii);
        bodies = dm(1:nMax);
        bds = sum(bodies~=0);
        if bds == 0
            continue
        end
        part2 = dm(nMax+1:bg);
        MC_POSITIONS(jj,1,ii) = dm(end);
        
        for kk = 1:bds
            MC_POSITIONS(jj,2+4*(kk-1):4*kk,ii) = dm(bg+1+3*(kk-1):bg+3*kk);
            MC_POSITIONS(jj,5+4*(kk-1),ii) = bodies(kk);
            shift = Site(1:3,Site(4,:)==bodies(kk));
            MC_POSITIONS(jj,2+4*(kk-1):4*kk,ii) = MC_POSITIONS(jj,2+4*(kk-1):4*kk,ii) - shift';
        end
        % Fill in remaining places with zeros
        MC_POSITIONS(jj,2+4*bds:end,ii) = 0;
        % Find the [ 0 0 0 S ] coordinate and move to front for eventual deletion
        part = zeros(nMax,4);
        for kk = 1:bds
            part(kk,:) = MC_POSITIONS(jj,2+4*(kk-1):1+4*kk,ii);
            if isequal(part(kk,1:3),[0 0 0]) && part(kk,4) == MC_POSITIONS(jj,1,ii)
                foundit = kk;
            end
        end
        dummy = MC_POSITIONS(jj,2:5,ii);
        MC_POSITIONS(jj,2:5,ii) = MC_POSITIONS(jj,2+4*(foundit-1):1+4*foundit,ii);
        MC_POSITIONS(jj,2+4*(foundit-1):1+4*foundit,ii) = dummy;
    end
end

% remove all unused space after the coordinates
MC_POSITIONS(:,2+4*nMax:end,:) =[];
% remove the central coordinates
MC_POSITIONS(:,2:5,:) = [];

%%%%%%%


%%% Write MC_POSITIONS to files
if exist('MC_POSITIONS','dir')==7
    rmdir 'MC_POSITIONS' s
    mkdir('MC_POSITIONS')
else
    mkdir('MC_POSITIONS')
end
cd './MC_POSITIONS'
for ii=1:sz3
    dlmwrite(num2str(ii),MC_POSITIONS(:,:,ii),'delimiter','\t','precision','%4.6g')
end
cd '../'

% Create a table from the sorted interactions with the heading created above
final_out = array2table(new_inter,'VariableNames',Headings);

fprintf('\n----------------------------------------------------------------\n\n')
fprintf('DONE!\nHere are the  n-body clusters that were found:\n\n')
disp(final_out)
fprintf('A total of %d unique clusters were found\n',num_interactions)
fprintf('These results have been written to "OUTPUT_INTERACTIONS.txt"\nRename as "INTERACTIONS.txt" in order to use as input to "COUNTS_GEN.m"\n')




%%% Write the output to file "OUTPUT_INTERACTIONS.txt"
fileID = fopen('OUTPUT_INTERACTIONS.txt','w');
fprintf(fileID,'   %s\t',Headings{1:end-1});
fprintf(fileID,'   %s\n',Headings{end});
fclose(fileID);
dlmwrite('OUTPUT_INTERACTIONS.txt',new_inter,'delimiter','\t','precision','%4.6g','-append')

fileID = fopen('CLUSTER_POSITIONS.txt','w');
fprintf(fileID,'   %s\t',Headings2{1:end-1});
fprintf(fileID,'   %s\n',Headings2{end});
fclose(fileID);
dlmwrite('CLUSTER_POSITIONS.txt',pos_new_inter,'delimiter','\t','precision','%4.6g','-append')

%Timing stuff
elapsed=toc;
inmin = elapsed/60;
fprintf('\nThis run took %9.2f seconds (or %3.2f min) to run.\n',elapsed,inmin)
fprintf('\n----------------------------------------------------------------\n')
fprintf('\n----------------------------------------------------------------\n\n')