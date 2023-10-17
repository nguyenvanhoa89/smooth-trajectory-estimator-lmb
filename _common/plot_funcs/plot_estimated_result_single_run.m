function plot_estimated_result_single_run(mc_results,mc_idx)
    model = mc_results.model;
    truth = mc_results.truth;
    
    est = mc_results.est{mc_idx};
    meas = mc_results.meas{mc_idx};

    mot_filter = strrep(model.p.Results.mot_filter,'_',' ');
    plot_truth(truth,1,truth.K); plot_tracks(est,1,length(est.X)); title(mot_filter);
    plot_card(truth,est); title(mot_filter);

%     plot_results(model,truth,meas,est,'plot_flag',true,'save_plot_flag',false);

end