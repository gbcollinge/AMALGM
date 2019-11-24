%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clean_configs.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


all_configs = dir('./configs');         % command 'dir' assigns each file a fileID starting from 3
num_configs = numel(all_configs)-2;     % since fileID 1 and 2 are STDIN and STDERR we want to discount these
confnum = 0;
for jj=1:num_configs                    % start working through the configurations...
    fprintf('CONFIGURATION #%d, %s:\n',jj, all_configs(jj+2).name)
    backhome = cd('./configs/');
    confnum = confnum + 1;
    analyzed_flag = 0;
    flag = 0;
    fID = fopen(all_configs(jj+2).name);
    tline = fgetl(fID);
    line=0;
    while ischar(tline)  
        if isempty(tline)
            tline = fgetl(fID);
            continue
        end       
        if strcmp(tline,'analyzed')                         % check if this configuration has already been analyzed
            flag = 1;       
            break
        else
        line = cell2mat(textscan(tline, '%f'))';
        dlmwrite('temp.txt',line,'precision','%12.10g','delimiter',' ','-append')
        tline = fgetl(fID);   
        end
    end   
    fclose(fID);
    if flag == 1    
        movefile('temp.txt',all_configs(jj+2).name)       
    end   
    if exist('temp.txt','file') ~=0
        delete('temp.txt')
    end
    cd(backhome)
end

