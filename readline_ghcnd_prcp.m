function [dayval]=readline_ghcnd_prcp(tline)

% info about line we're reading from http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% YEAR         12-15   Integer
% MONTH        16-17   Integer
% ELEMENT      18-21   Character
% VALUE1       22-26   Integer
% MFLAG1       27-27   Character
% QFLAG1       28-28   Character
% SFLAG1       29-29   Character
% VALUE2       30-34   Integer
% MFLAG2       35-35   Character
% QFLAG2       36-36   Character
% SFLAG2       37-37   Character
%   .           .          .
%   .           .          .
%   .           .          .
% VALUE31    262-266   Integer
% MFLAG31    267-267   Character
% QFLAG31    268-268   Character
% SFLAG31    269-269   Character
% ------------------------------

% example of line we're reading: 
%USW00093986191003PRCP    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6  127  6    0  6    0T 6   13  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0T 6    0  6    0  6


month=str2num(tline(16:17));
for day=1:31
    start=22+(day-1)*8;
% Changed trace method april 18 2018
    %    if strcmp(tline(start+4:start+5),'0T')
    %    dayval(day)=eps; % this is the trace value i'll use
    %else

    %disp(tline(start:start+4))
    %disp(tline(start+5:start+7))
        dayval(day)=str2num(tline(start:start+4));
        %end
    mflag(day)=tline(start+5);
    qflag(day)=tline(start+6);
    sflag(day)=tline(start+7);
end
dayval(dayval==-9999)=NaN;
dayval=dayval/10;

% get rid of anything that has a quality flag
if length(strfind(qflag,' '))<31
    days=1:31;
    days(strfind(qflag,' '))=[];
    %    disp(['a qualify flag! '  qflag(days)])
    dayval(days)=NaN;
    %                 break;break;break;break;
end

% check that 'Trace' is zero; then set to eps
trdays=strfind(mflag,'T');
if length(trdays)>0
    vtr=dayval(trdays);
    vtrnz=find(vtr~=0);
    if length(vtrnz)>0
        disp(['non-zero trace: ' num2str(vtr(vtrnz))])
    end
    %    dayval(trdays(vtr==0))=eps; % nnt
    dayval(trdays)=eps; % at
end

% Added april 18 2018

% check for 'missing presumed zero' days
mpzdays=strfind(mflag,'P');
if length(mpzdays)>0
    disp(['missing presumed zero: ' num2str(length(mpzdays)) ' days'])
    dayval(mpzdays)=0;
end

% check for source flag 'S' (use with caution) and omit
cdays=strfind(sflag,'S');
if length(cdays)>0
    %disp(['use with caution: ' num2str(length(cdays)) ' days'])
    %disp(tline)
    dayval(cdays)=NaN;
end

% check for cocorahs and omit
cdays=strfind(sflag,'N');
if length(cdays)>0
    disp(['cocorahs: ' num2str(length(cdays)) ' days'])
    dayval(cdays)=NaN;
end



