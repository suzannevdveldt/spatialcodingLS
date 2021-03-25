function msPlotFootprintsAcrossSessionsSV(ms, cell_registered_struct);

%% Alternatively:

%put here the first to last sessions you want to be looking at 

firstsession = 1;
lastsession = 6;
nsessions = length(firstsession:lastsession);
sessionID = [firstsession:lastsession]

  y = 1;
    n = 0;
    Y = 1;
    N = 0; 

cells_to_delete = [];

alligned_cells = cell_registered_struct.cell_to_index_map;
session_alligned_cells = alligned_cells(:,[firstsession:lastsession]);
session_alligned_cells = session_alligned_cells(all(session_alligned_cells,2),:);

for cell_id = 1:length(session_alligned_cells);   
            figure; 
            set(gcf,'Position',[100 100 1500 600])

            for session = 1:nsessions;
                if session == 1;

                subplot(1,nsessions,session);
                    cell_id_registered = session_alligned_cells(cell_id,session);
                    cell_outline =  permute(cell_registered_struct.spatial_footprints_corrected{sessionID(session),1}(cell_id_registered,:,:),[2 3 1]);
                    corrected_SFP = cell_registered_struct.spatial_footprints_corrected{sessionID(session),1};
                    h = max(permute(corrected_SFP,[2 3 1]),[],3);
                    cell_outline(cell_outline > 0) = 1.5*(max(h(:)));
                    outlined_SFPs = h + cell_outline;
                    imagesc(outlined_SFPs);
                    colormap(gca,'jet')
                    hold on
                       
                else 
                subplot(1,nsessions,session);
                cell_id_registered = session_alligned_cells(cell_id,session);
                
                    subplot(1,nsessions,session);
                    cell_outline =  permute(cell_registered_struct.spatial_footprints_corrected{sessionID(session),1}(cell_id_registered,:,:),[2 3 1]);
                    corrected_SFP = cell_registered_struct.spatial_footprints_corrected{sessionID(session),1};
                    h = max(permute(corrected_SFP,[2 3 1]),[],3);
                    cell_outline(cell_outline > 0) = 1.5*(max(h(:)));

                    outlined_SFPs = h + cell_outline;
                    imagesc(outlined_SFPs);
                    colormap(gca,'jet')

                end
            end   % to go to next session
     
                        cell_id
                        prompt = 'Keep this cell? type Y for yes and N for no        ';
                        x = input(prompt);
                            if x == 0;
                            cells_to_delete(end+1) = cell_id; 
                            else
                            end

                        %pause
                        close all 
   
end

    session_alligned_cells([cells_to_delete],:) = []; 
   % alligned_cells = alligned_cells(:,[1 3]); 
    
save('session_alligned_cells.mat', 'session_alligned_cells');
end 
 
