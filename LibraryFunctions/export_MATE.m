function [rup] = export_MATE(rup)

% By Victor Hernández (victorh@hi.is). July 2025
% Function to export the kinematic source model in SPEED format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = rup.glo_coor(:,1);
    y = rup.glo_coor(:,2);
    z = rup.glo_coor(:,3);
    slip = rup.slip_deterministic(:);
    tau = rup.risT.dist{1,1}(:);
    trup = rup.rupT.dist{1,1}(:);
    pt = rup.PT.dist{1,1}(:);
    SeisMom = rup.Mo_vec;

    npts = length(SeisMom);
    
    % Writing .mate file
fid = fopen(fullfile([cd,'\',rup.Filename,'\',rup.Filename,'.mate']),'w');

isism = 0;

if strcmp(rup.svf.type,'exp') %Different format depending on the type of STF
    for ipt = 1:npts
        if slip(ipt)>0
            isism = isism + 1;
            SV=rup.slip_vectors(ipt,:);
            NV=rup.normal_vectors(ipt,:);
 
            fprintf(fid,'SISM  %d  %d  %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e     %+13.7e    %+13.7e    %+13.7e \n',...
            10, 2, rup.hypo_coor(1), rup.hypo_coor(2), rup.hypo_coor(3), x(ipt), y(ipt), z(ipt), SV(1), SV(2), SV(3), NV(1), NV(2), NV(3), trup(ipt), SeisMom(ipt), tau(ipt));
        end
    end
else
        for ipt = 1:npts
        if slip(ipt)>0
            isism = isism + 1;
            SV=rup.slip_vectors(ipt,:);
            NV=rup.normal_vectors(ipt,:);
 
            fprintf(fid,'SISM  %d  %d  %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e     %+13.7e    %+13.7e    %+13.7e  %+13.7e \n',...
            10, 2, rup.hypo_coor(1), rup.hypo_coor(2), rup.hypo_coor(3), x(ipt), y(ipt), z(ipt), SV(1), SV(2), SV(3), NV(1), NV(2), NV(3), trup(ipt), SeisMom(ipt), tau(ipt), pt(ipt));
        end
        end
end
    fclose(fid);

end