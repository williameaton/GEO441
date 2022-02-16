clc
clear all 
close all 

%axis equal
a=zeros(500,400);
fid= fopen('2homorough.out','r');
%cl=[-4.0:0.2:4.0];
for i=1:100;
   a(:,:)=fread(fid,[500,400],'float32');
   pcolor(a(:,:)); shading flat;
   colorbar('horiz');
   caxis([-0.5 0.5]);
	
	
   pause(0.1);
end
fclose(fid);
