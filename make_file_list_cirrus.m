% Script for generating a list of .img files contained within a directory.

% ** Type directory to search here **
base_path = 'C:/OCTData/Cirrus';
save_filename = 'file_list_cirrus.txt';


%% code to get list

root = fullfile([base_path,'/**/*Macular Cube*_z.img']);
files = rdir(root);

fid = fopen(save_filename,'w');
fprintf(fid,'%s',files(1).name);
for i = 2:length(files)
    fprintf(fid,'\n%s',files(i).name);
end
fclose(fid);