cromshell submit CollectPEMetricsForCPX.wdl CollectPEMetrics.CPX_under_1perc.json ../monitor_options.gnomad.json dependencies.zip

cromshell submit CollectPEMetricsForCPX.wdl CollectPEMetrics.CPX_1to10perc.json ../monitor_options.gnomad.json dependencies.zip 

cromshell submit CollectPEMetricsForCPX.wdl CollectPEMetrics.CPX_over10perc.json ../monitor_options.gnomad.json dependencies.zip

cromshell submit ReviseVcf.wdl example_json.ReviseVcf.gnomAD_V3.json ../monitor_options.gnomad.json dependencies.zip 

cromshell submit ReviseVcf.wdl example_json.ReviseVcf.gnomAD_V3.chr2.json ../monitor_options.gnomad.json dependencies.zip