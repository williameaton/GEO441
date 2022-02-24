clc
clear all 
close all 

z = zeros(1,101);    % maintain height of 0 m
x = 0:1000:100000;   % 100 km in 1 km increments
y = x;  
origin = [42.3648, -71.0214, 10.0];



local2latlon(x,y,z,origin)

%axis equal
a=zeros(500,400);
fid= fopen('2heterorough.out','r');
%cl=[-4.0:0.2:4.0];
for i=1:100;
   a(:,:)=fread(fid,[500,400],'float32');
   pcolor(a(:,:)); shading flat;
   colorbar('horiz');
   caxis([-1 1]);
	
	
   pause(0.1);
end
fclose(fid);
