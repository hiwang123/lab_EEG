
function y=filter_idealfilter(raw,lbound,ubound,fz)

len=size(raw(:));
ts=timeseries(raw',(1:len)/fz);     % data -> timeseries class
fl = idealfilter(ts, [lbound ubound], 'pass');    % pass:keep []  notch:remove[]

y=fl.Data;