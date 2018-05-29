function V = queryPP(name,parameter);

%name='c:/img/e1720/data/gems_01.img'
%parameter='name'
%----------------------------------------
%function queryPP
%finds a parameter value in procpar
%----------------------------------------
%Usage V = queryPP(name,par);
%
%Input:
%name = name of FID directory without the .fid extension
%par  = name of parameter
%
%Output:
%V    = value of parameter - possibly a vector or matrix (strings)
%
%Examples:
%V = queryPP('sems_ms','gro');
%
%
%----------------------------------------
% Maj Hedehus, Varian, Inc., Oct 2001.
%----------------------------------------


if exist('name') == 0
    name = input('DIR name : ','s');
end

if exist('parameter') == 0
    parameter = input('Which parameter: ','s');
end


fullname = sprintf('%s%cprocpar',name,'/'); %something funky with the backslash going on...

fid = fopen(fullname,'r');
if fid == -1
    str = sprintf('Can not open file %s',fullname);
    error(str);
end

par  = fscanf(fid,'%s',1);

while (~strcmp(par,parameter)) & ~feof(fid)
    type  = fscanf(fid,'%d',1);
    fgetl(fid); %skip rest of line
    nvals = fscanf(fid,'%d',1);
    fgetl(fid); %skip rest of line
    if type == 2 %string, the remaining values are on seperate lines
        for n = 1:nvals-1
            fgetl(fid);
        end
    end
    fgetl(fid);  %skip "line 3"
    
    par = fscanf(fid,'%s',1);
end

if feof(fid)
    V = 'parameter not found';
    return
end


type  = fscanf(fid,'%d',1);
fgetl(fid); %skip rest of line
nvals = fscanf(fid,'%d',1);

switch type
case 1  % float
    V = fscanf(fid,'%f',nvals);
case 2  % string
    L = fgetl(fid);
    if (strcmp(parameter,'name') || strcmp(parameter,'comment'))
        sm=L(2:end);
    else
        sm = sscanf(L,'%s',1); %for some odd reason, fscanf(fid,..) doesn't work here
    end
    for n = 1:nvals-1
        L = fgetl(fid);
        sm = sscanf(L,'%s',1);
        sm = str2mat(sm); %,s);
        V{n} = sm(2:length(sm)-1); %skip " " quotes
    end        
    if (nvals==1)
        V=sm(2:length(sm)-1);
    end
    
case 3  % delay
    V = fscanf(fid,'%f',nvals);
case 4  % flag
    L = fgetl(fid);
    sm = sscanf(L,'%s',1); %for some odd reason, fscanf(fid,..) doesn't work here
    for n = 1:nvals-1
        L = fgetl(fid);
        sm = sscanf(L,'%s',1);
        sm = str2mat(sm)
    end        
    V = sm;
case 5  % frequency
    V = fscanf(fid,'%f',nvals);
case 6  % pulse
    V = fscanf(fid,'%f',nvals)*1e6;
case 7  % integer
    V = fscanf(fid,'%d',nvals);
otherwise
    V = 'not defined';
end

fclose(fid);