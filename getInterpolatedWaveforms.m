function varargout=getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName)
%
% waveforms= getInterpolatedWaveforms(varLabels,category)
%
% Generates interpolated waveforms from unevenly sampled data 
% in category.
%
% varLabeLs Nx1 cell array of the physiological labels to category of the waveforms to be
% interpolated
%
% category Mx3 cell array 
%
% Ts 1x1 sampling interval (in hours) from which to interpolate data
%
% outVarName Nx1 cell array for the names of the new interpolated data
%
% Output is multiple variables, one for each element of varLabels see 
% example below:
%
% %Example 1:
% outVarName={'lact','map','hr','urine'};
%[lact,map,hr,urine]=getInterpolatedWaveforms(varLabels,category,Ts,outVarName);

NvarName=length(varLabels);

nb=round(0.5/Ts); %Filter waveforms with half an hour moving average
b=ones(nb,1)./nb;
for n=1:NvarName
    
    ind=strcmp(category,varLabels{n});
    x=[tm(ind) val(ind)];
    if(length(x)<3)
        eval([outVarName{n} '=[];'])
        continue;
    end
    x=sortrows(x,1);
    del=find(isnan(x(:,1))==1);
    x(del,:)=[];
    del=find(x(:,2)==0);
    if(~isempty(del))
        x(del,:)=[];
    end
    %When dealing with urine, remove first 2 values because they maybe
    %cumulative over the operating room stay
    if(strcmp(outVarName{n},'urine') && length(x(:,1)))
            x(1:2,:)=[];
    end
    
    if(length(x)<3)
        %Not enough minimum datapoints, set waveform to empty
        eval([outVarName{n} '=[];'])
        continue;
    end
    
    %Only keep the longest segment that has a minimum of P points
    %per hour for TmDuration hours
    isGood=hasFrequentHourlyPoints(x,P,TmDuration);
    
    x=hourly_median(x);
    
    %Interpolate waveforms
    y=[x(1,1):Ts:x(end,1)]';
    y(:,2)=interp1(x(:,1),x(:,2),y,'linear');
    
    %Filter the waveforms through a moving average
    y(:,2)=filtfilt(b,1,y(:,2));
    
    %Set y to the time series being analyzed
    eval([outVarName{n} '=y;'])
end
for n=1:nargout
    eval(['varargout{n}=' outVarName{n} ';'])
end
