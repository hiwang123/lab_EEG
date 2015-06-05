function [TD LE]=cutting(td,ST,ET,sr)
    le=length(td);
    st=ST*sr;
    et=ET*sr;
    if isempty(ET)
        et = le;
    end
    if et>le 
        et=le;
    end
    if st<=0
        st=1;
    end
    TD=td(:,st:et-1);
    LE=(1:length(TD))/sr;
end