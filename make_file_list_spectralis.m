% Script for generating a list of .vol files contained within a directory.

% ** Type directory to search here **
base_path = 'C:/OCTData/Spectralis';
save_filename = 'file_list_spectralis.txt';


%% code to get list

root = fullfile([base_path,'/**/*.vol']);
files = rdir(root);

fid = fopen(save_filename,'w');
fprintf(fid,'%s',files(1).name);
for i = 2:length(files)
    fprintf(fid,'\n%s',files(i).name);
end
fclose(fid);