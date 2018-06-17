function [ AC ] = parseRunCaseFile( filename )
% Copyright 2012 Joseph Moster
% Adapted to parse all necessary data by RJ Gritter, September 2012
% -fileLen mod done
% Adapted again by Ben Brelje 2012
% Requires the aero design toolbox folder

file = textread(filename, '%s', 'delimiter', '\n','whitespace', '');
i=1;

fileLen = length(file);
while i<fileLen
    str = char(file(i));
    header = regexpi(str, 'Run case');
    
    if(~isempty(header))
        
        %Skip remainder of header
        i=i+1;
        while(~isempty(char(file(i))))
            i=i+1;
        end
        
        %Read off values
        AC.Alpha   =  findValue(file,'Alpha', [i,fileLen]);
        AC.CLtot   =  findValue(file,'CLtot', [i,fileLen]);
        AC.CDtot   =  findValue(file,'CDtot', [i,fileLen]);
        AC.e       =  findValue(file,'e =', [i,fileLen]);

    end
    i=i+1;
end
end