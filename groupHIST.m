function [grpPDF,grp_elORDER,grpHIST]=groupHIST(GRP,binCEN)
% GROUPHIST(GRP,binCEN) bins the group data, returns the distribution of
% group data in each channel. the bin center is specified as a vector in
% binCEN. Resulting histogram each has bin number as the length of binCEN.
%   
%   <INPUT> 
%      GRP: group data as a m-channel-by-n-time-point matrix. 
%   binCEN: specified itself as a vector of bin centers, the length of which
%           is the number of bins you would like to bin your data.
% 
%   <OUTPUT>
%        grpPDF: 
%   grp_elORDER:
%       grpHIST:
% 
%   <EXAMPLE>
%           00 binCEN = linspace(<MIN>,<MAX>,<# of bins>);
%           01 [<outputs>] = groupHIST(<data>,binCEN); 
%
% 
% *** THIS SUBROUTINE IS IMPLEMENTED AS A PART OF `TFCMI' ***
% 
% Yen Yu <higgs.bose@gmail.com>, IBRL, VGH Taipei
% 2007/02/13
% rev. #1 2007/02/13
% original version from `jointp.m' by Chan-Chuan Chen, 2003.04


[nChannel,nelTimeCourse]=size( GRP );
nBins = length(binCEN);
grpHIST = hist( GRP',binCEN );
grpHIST = grpHIST';
grpHIST( find(grpHIST==0) ) = 1e-16;
grpPDF  = grpHIST/nelTimeCourse;
grp_elORDER = zeros(nChannel,nelTimeCourse);

for CH=1:nChannel
    % Determine the order of elements in each channel with respect to the
    % histogram order
	for eTC=1:nelTimeCourse
    		ordTmp=abs( binCEN - GRP(CH,eTC) );
    		[dummy, iTmp]=min( ordTmp );
    		grp_elORDER(CH,eTC)=iTmp;
	end

end
