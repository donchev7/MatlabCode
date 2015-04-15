function filewrite(myfile, mytext)
% Writes a basic text file to some file name analogous to the matlab 
% "fileread" function.
% Usage:
% filewrite(filename, 'the dog ate the fish');

% Jade Q. Wang
% http://www.jadism.com
% 07.07.10

    fid = fopen(myfile,'w');
    fprintf(fid, '%s', mytext);
    fclose(fid);