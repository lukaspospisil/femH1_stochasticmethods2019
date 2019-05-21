function plot_eeg_movie( P, M, X, T_to_plot )

% create a mesh grid on which we will plot
[XI,YI]=meshgrid([-8:0.1:8],[-8:0.1:6]);

fig_movie = figure;
col=[min(min(X)) max(max(X))];

disp(['- movie of the first ' num2str(length(T_to_plot)) ' instances off EEG measurement'])
disp( '   [movie is paused, press a key to process to next frame]')

%% show the 2D-movie of the first "T" time instances of the 
%% EEG measurement 
for tt=1:length(T_to_plot)
   t = T_to_plot(tt);

   clf
   ZI(:,:,t) = griddata([P(:,1);M(:,1)],...
            [P(:,2);M(:,2)],[X(:,t);zeros(size(M,1),1)],XI,YI);
   contourf(XI,YI,ZI(:,:,t));caxis(col);
   title(['    t = ' num2str(t)])
   disp(['    (number of frame: ' num2str(t) ' / ' num2str(length(T_to_plot)) ')'])
   pause
 end

close(fig_movie)


end

