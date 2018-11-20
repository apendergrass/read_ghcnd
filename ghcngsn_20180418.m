diri='/datalocal/ccpd/apgrass/ghcn/2018-04-18/'

ghcninv=importdata([diri 'ghcnd-inventory.txt']);
%         data: [nstationx2 double]
%     textdata: {nstationx4 cell}

% FORMAT OF "ghcnd-inventory.txt"
% 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEMENT      32-35   Character
% FIRSTYEAR    37-40   Integer
% LASTYEAR     42-45   Integer
% ------------------------------

% find stations with precip data ending in at least 2014 with at least 50 years of data
%goodones=find((ghcninv.data(:,2)>2014)&(diff(ghcninv.data,1,2)>49)&strcmp('PRCP',ghcninv.textdata(:,4)));
%goodones=find((ghcninv.data(:,2)>2014)&(diff(ghcninv.data,1,2)>49)&strcmp('PRCP',ghcninv.textdata(:,4))&(~strncmp('US1',ghcninv.textdata(:,1),3)));
%~strncmp('US1',ghcninv.textdata(:,1),3) %cocorahs stations


%goodones=find((ghcninv.data(:,2)>2014)&(diff(ghcninv.data,1,2)>49)&strcmp('PRCP',ghcninv.textdata(:,4))&(~strncmp('US1',ghcninv.textdata(:,1),3)));
%&(diff(ghcninv.data,1,2)>49)
%goodones=find(strcmp('PRCP',ghcninv.textdata(:,4))&(~strncmp('US1',ghcninv.textdata(:,1),3))&(diff(ghcninv.data,1,2)>19));  %
goodones=find(strcmp('PRCP',ghcninv.textdata(:,4))&(~strncmp('US1',ghcninv.textdata(:,1),3)));
% precip data, not from cocorahs
                                                      % excluding
                                                      % cocorahs
                                                      % 966 with no
                                                      % time restrictions
years=ghcninv.data(goodones,:);

% strcmp('PRCP',ghcninv.textdata(:,4))

figure(1);clf
invids=ghcninv.textdata(goodones,1);
statlats=ghcninv.textdata(goodones,2);
statlons=ghcninv.textdata(goodones,3);

for i=1:length(goodones)
    lats(i)=str2num(statlats{i});
    lons(i)=str2num(statlons{i});
end

hold on 
plot(lons,lats,' kx')
plot(lons,lats,' rx')


dateres=length(invids) % all prcp: 104270 
                       % excluding cocorahs: 72627
                       % with date restrictions: 12111

ghcnstations=importdata([diri 'ghcnd-stations.txt']);

% FORMAT OF "ghcnd-stations.txt" 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEVATION    32-37   Real
% STATE        39-40   Character
% NAME         42-71   Character
% GSN FLAG     73-75   Character
% HCN/CRN FLAG 77-79   Character
% WMO ID       81-85   Character
% ------------------------------

for i=1:size(ghcnstations.textdata,1)
    t=textscan(ghcnstations.textdata{i},'%s %f %f %f %31c %s %s');
    sid{i}=t{1};
    elev{i}=t{4};
    elev{i}=t{4};
    name{i}=t{5};
    gsn{i}=t{6};
    theend{i}=t{7};
end

for i=1:size(ghcnstations.textdata,1)
    t=textscan(ghcnstations.textdata{i},'%s %f %f %f %31c %s %s');
end



count=0;
for i=1:length(theend)
    if ~isempty(theend{i}); 
        count=count+1;
        notempty(count)=i;
    end
end

% just checking.
% count=0;
% for i=1:length(theend)
%     if ~isempty(strfind(name{i},'GSN')); 
%         count=count+1;
%         notempty(count)=i;
%     end
% end
% 
% count=0;
% for i=1:length(theend)
%     if ~isempty(strfind(name{i},'GSN')); 
%         count=count+1;
%         notempty(count)=i;
%     end
% end

count=0;
for i=1:length(theend)
    if strcmp('GSN',gsn{i});
        count=count+1;
        isgsn(count)=i;
    end
end
isgsn=isgsn(:);

for i=1:length(notempty)
    if strcmp('GSN',theend{notempty(i)});
        isgsn(end+1)=notempty(i);
    end
end

isgsn=sort(isgsn);

sid=sid(isgsn);
elev=elev(isgsn);

count=0;
for i=1:length(isgsn)    
    in=find(strcmp(sid{i},invids));
    if ~isempty(in)
        count=count+1;
        all2gsn(count)=in;
    end
end

gsnlats=lats(all2gsn);
gsnlons=lons(all2gsn);
gsnids=invids(all2gsn);

plot(lons,lats,' kx',gsnlons,gsnlats,' rx')
legend('all ghcnd','gsn','location','northoutside')
%xlim([-140 -60])
%ylim([20 60])

xlim([-180 180])
ylim([-90 90])
set(gcf,'units','centimeters','paperpositionmode','auto');
set(gcf,'position',[10 14 8.4 8.4]); %% one col: 8.4 cm. two col 16.9.  max height 23.7.

%epswrite('ghcndgsnmap.eps')
%epsfixfonts('ghcndgsnmap.eps')

%119 poc% ls -1 > fnames.txt
%119 poc% pwd
% 123 poc% ls -1 | wc -l

% check to make sure we really have all these files 
fnames=importdata([diri '/gsn/ghcnd_gsn/fnames.txt']);

for i=1:(size(fnames,1)-1)
    t=textscan(fnames{i},'%11c %s');
    fileids{i}=t{1};
end

count=0;
missingfile=[]
for i=1:length(gsnids)
    if isempty(find(strcmp(gsnids(i),fileids(:))))
        count=count+1;
        missingfile(count)=i;
    end
end
missingfile


gsnids(missingfile)=[];
gsnlats(missingfile)=[];
gsnlons(missingfile)=[];
all2gsn(missingfile)=[];

% sort(years(all2gsn,1))
years=ghcninv.data(goodones,:);
gsnstartyears=years(all2gsn,1);
gsnendyears=years(all2gsn,2);
firststartyear=min(gsnstartyears)
laststartyear=max(gsnstartyears)

startyears=firststartyear:laststartyear;
for i=1:length(startyears)
    nstaty(i)=sum(gsnstartyears<=startyears(i));
end

figure(2);
plot(startyears,nstaty,'ok-')


startyear=1893; % this year there's a big jump in the number of
                % stations online month=str2num(tline(16:17));
endyear=str2num(datestr(now,10))-1;

% 	   The five core elements are:
% 
%            PRCP = Precipitation (tenths of mm)
%    	   SNOW = Snowfall (mm)
% 	   SNWD = Snow depth (mm)
%            TMAX = Maximum temperature (tenths of degrees C)
%            TMIN = Minimum temperature (tenths of degrees C)

prcp=NaN([length(gsnids) length(startyear:endyear+1) 12 31]);
%snow=NaN([length(gsnids) length(startyear:endyear) 12 31]);
%snwd=NaN([length(gsnids) length(startyear:endyear) 12 31]);
%tmax=NaN([length(gsnids) length(startyear:endyear) 12 31]);
%tmin=NaN([length(gsnids) length(startyear:endyear) 12 31]);

for stati=1:length(gsnids)
    disp([num2str(stati) ' ' gsnids{stati}])
    fid=fopen(strcat(diri,'gsn/ghcnd_gsn/',gsnids{stati},'.dly'));
    tline=fgetl(fid);
    while ischar(tline)
        year=str2num(tline(12:15));
        if year>(startyear-1)
            month=str2num(tline(16:17));
            variable=tline(18:21);
            switch variable
                case 'PRCP' %            if strcmp('PRCP',tline(18:21))
                    prcp(stati,year-(startyear-1),month,:)=readline_ghcnd_prcp(tline);
                    % readline_ghcnd_prcp_nnt does the same thing
                    % but treats traces on days with something
                    % other than zero precip as whatever the
                    % reading was (it's usually missing)
            end
        end
        tline=fgetl(fid);
    end
    fclose(fid);
end

prcpraw=prcp; 


pnz=prcp;
pnz(prcp>0)=1;
prcp(prcp==eps)=0;



years=startyear:endyear+1;

save gsndailydata_20180418.mat prcp pnz gsnids gsnlats gsnlons ...
    years gsnstartyears gsnendyears

                                                                      


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % write a netcdf file of this stuff.




% filename=[diri 'ghcnd_gsn_prcp.nc'];                                                                                                   


% prcp(find((prcp==0)&(pnz==1)))=1e-6;

% %%% write years
% nccreate(filename,'year','Dimensions', {'year', length(years)}); %???
% ncwrite(filename,'year',years);                                                                                                                        

%     ncwriteatt(filename,'/','creation_date',datestr(now));
%     ncwriteatt(filename,'/','description','GHCN Daily Historical Climate Network core variables, 1893-present. Any flagged data is omitted. http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/');

% %%% write station IDs
% nccreate(filename,'ID','Dimensions', {'cols',11,'ID',length(gsnids)},'DataType','char');                                                   
% ncwrite(filename,'ID',strvcat(gsnids)'); %???

% %%% write station lat/lons
% nccreate(filename,'lat','Dimensions', {'station', length(gsnids)});                                                                                
% nccreate(filename,'lon','Dimensions', {'station', length(gsnids)});                                                                                
% % nccreate(filename,'year','Dimensions', {'year', length(years)});                                                                                

% ncwrite(filename,'lat',gsnlats);                                                                                                                        
% ncwrite(filename,'lon',gsnlons);                                                                                                                        

% %%% write data 
% nccreate(filename,'prcp','Dimensions', {'station', length(gsnids),'year',length(years),'Month',12,'Day',31});                                                                                

% ncwriteatt(filename,'prcp','longname','Precipitation');
% ncwriteatt(filename,'prcp','units','mm');
% ncwriteatt(filename,'prcp','description','Precipitation. Trace set to 1e-6');

% ncwrite(filename,'prcp',prcp);                                                                                                                        


