function jPDF=jpdf4CMI(xORDER,yORDER,nBINs)
% JPDF4CMI(XORDER,YORDER,NBINS)  calculates the joint probability for two
% histogram for the calculation of cross mutual information (cMI). 
% 
%   <DEPENDS> groupHIST.m
% 
%   <INPUT>
%       xORDER,yORDER: element orders from groupHIST.m
%               nBINs: number of bins used to bin to generate the
%                      probability density function.
% 
%   <OUTPUT>
%       jPDF: 
% 
%   <EXAMPLE>
% 
% 
% 
%  
% 
% 
% 
% *** THIS SUBROUTINE IS IMPLEMENTED AS A PART OF `TFCMI' ***
% 
% Yen Yu <higgs.bose@gmail.com>, IBRL, VGH Taipei
% 2007/02/13
% rev. #1 2007/02/13
% original version from `jointp3.m' by Chan-Chuan Chen, 2003.04



nelTimeCourse=size(xORDER,2);
% jPDF = ones(nBINs,nBINs)*1e-16; % FIXME: to be deleted
jPDF = zeros(nBINs,nBINs);
for t=1:nelTimeCourse
    % counts and accumulates where events coincide.
    jPDF( xORDER(t),yORDER(t) )= jPDF( xORDER(t),yORDER(t) )+1;
    
end

jPDF = jPDF/nelTimeCourse;
jPDF( find(jPDF==0 ) ) = 1e-16;