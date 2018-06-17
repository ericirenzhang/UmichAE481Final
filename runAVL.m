function [AVLdata] = runAVL()

%% File Locations
avl_loc = '.\AVL\avl.exe';
filename = '.\AVL\airplane.avl';
runname = '.\AVL\runavl';
dataname = '.\AVL\DATA\\TEST_JM';
casefile = '.\AVL\case';

%% Write to run case file


%% Create file to execute AVL
[status, result] = system(strcat('del ',runname, '.run'));
fid = fopen(strcat(runname,'.run'),'W');
fprintf(fid,'LOAD %s\n',filename);
fprintf(fid, '%s\n', 'OPER');
fprintf(fid, '%s\n','F');
fprintf(fid, '%s\n',casefile);
fprintf(fid, '%s\n','X');
fprintf(fid, '%s\n','FT');
fprintf(fid, '%s\n',dataname);
fprintf(fid, '%s\n','O');
fprintf(fid, '%s\n','');
fprintf(fid, '%s\n','');
fprintf(fid, 'Quit\n');
fclose(fid);

%% Execute AVL
fid = system(strcat(avl_loc,' < ',runname,'.run>null'));
% fprintf('\n');

AVLdata = parseRunCaseFile(dataname);