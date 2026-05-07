%Simfldr_path = '..\Examples\Run_Sims_a0006\expected_outputs\a-0006-ap=1.614011963mm';
Simfldr_path = '..\Examples\Run_Sims_p0004\expected_outputs\p-0004-ap=0.700283325mm';

Q_to_display = 'ext';
figs_path = '';
plot_cells = 1;
plot_MIDAS_Qid_fromSimPath(Simfldr_path,Q_to_display,figs_path,plot_cells);
Q_to_display = 'sca';
plot_cells = 0;
plot_MIDAS_Qid_fromSimPath(Simfldr_path,Q_to_display,figs_path,plot_cells);
Q_to_display = 'bks';
plot_MIDAS_Qid_fromSimPath(Simfldr_path,Q_to_display,figs_path,plot_cells);
