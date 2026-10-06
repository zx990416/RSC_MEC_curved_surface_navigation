function v2 = smooth_data_sep_nan(v, smoothwidth, smooth_min_bin)
if(~exist('smoothwidth','var'))
    smoothwidth = 5;  % 如果未出现该变量，则对其进行赋值
end
if(~exist('smooth_min_bin','var'))
    smooth_min_bin = 4;  % 如果未出现该变量，则对其进行赋值
end
v2 = nan(length(v),1);
M = [v,(1:length(v))'];
idx = all(isnan(v),2);
idx = isnan(v);
idy = 1+cumsum(idx);
idz = 1:size(v,1);
C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});
C = C(~cellfun('isempty', C));
ncell = length(C);
for icell = 1:ncell
    temp = C{icell};
    if size(temp,1)<smoothwidth
        continue
    end
    tempv = smooth(temp(:,1), smoothwidth);
    v2(temp(:,2)) = tempv;
end

end
