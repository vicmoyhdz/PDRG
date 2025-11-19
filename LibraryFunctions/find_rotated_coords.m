function [x_global,y_global] = find_rotated_coords(x_local,y_local,Theta)

% Srihari Sangaraju - 11-Nov-2018

% Used to find the Global coordinates of (x_local,y_local).
% where (x_local,y_local) are with respect to a local coordinate system with
% Origin as (0,0) and local_x_axis makes an angle "Theta" with respect to
% global_x_axis in anticlockwise direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_global = zeros(size(x_local));
y_global = x_global;

 for i=1:length(x_local(:))
    Radi = sqrt((x_local(i)^2)+(y_local(i)^2));
    rto = y_local(i)/x_local(i);

    if sign(rto) == 1 && x_local(i)>0         % 1st Quadrant in Local coordinates
        phi = atand(rto);
    elseif sign(rto) == -1 && x_local(i)<0    % 2nd Quadrant in Local coordinates
        phi = 180 + atand(rto);
    elseif sign(rto) == 1 && x_local(i)<0     % 3rd Quadrant in Local coordinates
        phi = 180 + atand(rto);
    elseif sign(rto) == -1 && x_local(i)>0    % 4th Quadrant in Local coordinates
        phi = atand(rto);
    elseif sign(rto) == 0 && x_local(i)>0
        phi = 0;
    elseif sign(rto) == 0 && x_local(i)<0
        phi = 180;
    elseif sign(rto) == -1 && x_local(i)==0
        phi = 270;
    elseif sign(rto) == 1 && x_local(i)==0
        phi = 90;
    elseif Radi ==0
        phi = 0;
    end

    x_global(i) = Radi*cosd(phi+ Theta);
    y_global(i) = Radi*sind(phi+ Theta);
 end
end

