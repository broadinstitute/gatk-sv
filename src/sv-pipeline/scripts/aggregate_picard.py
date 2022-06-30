def read_in(file_name, label):
    started = False
    newlist = []
    wgsfile = file_name
    
    with tfio.gfile.GFile(file_name, "r") as inp:
        for line in inp:
            intable = line.startswith(label)
            if intable == True:
                started = True
                continue
            if started:
                if line != '\n':
                    newlist.append(line.rstrip('\n'))
            if started and line == '\n':
                return newlist
                break
                
def newtable(file_name, label):
    columnnames = read_in(file_name, label)[0].split("\t")
    values = read_in(file_name, label)[1].split("\t")
    dataframe = pd.DataFrame(columnnames).T
    dataframe.columns = columnnames
    dataframe_length = len(dataframe)
    dataframe.loc[dataframe_length] = values
    dataframe = dataframe.drop(0)
    return dataframe

def addrows(file_name, label):
    onerow = newtable(file_name, label)
    for i in range(len(read_in(file_name, label))):
        nextrow = read_in(file_name, label)[i].split("\t")
        onerow.loc[i] = nextrow
    onerow = onerow.drop(0)
    return onerow

def concatrows(table, suffix):
    newdata = pd.DataFrame()
    for row in range(1, len(table) + 1):
        addrow = pd.DataFrame(table.loc[row]).T.reset_index(drop=True)
        addrow = addrow.add_suffix(suffix[row-1])
        newdata = pd.concat([newdata, addrow], axis = 1)
    return newdata

def window_metrics_mean(table, window):
    dataframe1 = pd.DataFrame(table[0:window]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    dataframe2 = pd.DataFrame(table[window:window*2]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    dataframe3 = pd.DataFrame(table[window*2:window*3]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    dataframe4 = pd.DataFrame(table[window*3:window*4]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    dataframe5 = pd.DataFrame(table[window*4:window*5]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    dataframe6 = pd.DataFrame(table[window*5:window*6]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    
    bddf1mean = dataframe1.mean()
    bddf2mean = dataframe2.mean()
    bddf3mean = dataframe3.mean()
    bddf4mean = dataframe4.mean()
    bddf5mean = dataframe5.mean()
    bddf6mean = dataframe6.mean()
    
    bddf1meanvalues = pd.DataFrame(bddf1mean).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window))
    bddf2meanvalues = pd.DataFrame(bddf2mean).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*2))
    bddf3meanvalues = pd.DataFrame(bddf3mean).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*3))
    bddf4meanvalues = pd.DataFrame(bddf4mean).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*4))
    bddf5meanvalues = pd.DataFrame(bddf5mean).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*5))
    bddf6meanvalues = pd.DataFrame(bddf6mean).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*6))
    
    BDMresult = pd.concat([bddf1meanvalues, bddf2meanvalues, bddf3meanvalues, bddf4meanvalues, bddf5meanvalues, bddf6meanvalues], axis = 1)
    return BDMresult

def window_metrics_std(table, window):
    mqdf1 = pd.DataFrame(table[0:window]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    mqdf2 = pd.DataFrame(table[window:window*2]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    mqdf3 = pd.DataFrame(table[window*2:window*3]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    mqdf4 = pd.DataFrame(table[window*3:window*4]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    mqdf5 = pd.DataFrame(table[window*4:window*5]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    mqdf6 = pd.DataFrame(table[window*5:window*6]).reset_index(drop = True).apply(pd.to_numeric, errors='ignore')
    
    mqdf1std = mqdf1.std().rename({'MEAN_QUALITY':'SD'})
    mqdf1stdtable = pd.DataFrame(mqdf1std).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window))
    mqdf2std = mqdf2.std().rename({'MEAN_QUALITY':'SD'})
    mqdf2stdtable = pd.DataFrame(mqdf2std).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*2))
    mqdf3std = mqdf3.std().rename({'MEAN_QUALITY':'SD'})
    mqdf3stdtable = pd.DataFrame(mqdf3std).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*3))
    mqdf4std = mqdf4.std().rename({'MEAN_QUALITY':'SD'})
    mqdf4stdtable = pd.DataFrame(mqdf4std).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*4))
    mqdf5std = mqdf5.std().rename({'MEAN_QUALITY':'SD'})
    mqdf5stdtable = pd.DataFrame(mqdf5std).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*5))
    mqdf6std = mqdf6.std().rename({'MEAN_QUALITY':'SD'})
    mqdf6stdtable = pd.DataFrame(mqdf6std).T.drop('CYCLE', axis = 1).add_suffix('_'+str(window*6))

    MQMresult = pd.concat([mqdf1stdtable, mqdf2stdtable, mqdf3stdtable, mqdf4stdtable, mqdf5stdtable, mqdf6stdtable], axis = 1)
    return MQMresult

def sequencing_artifact_metrics(table):
    sadfone = pd.DataFrame(table.loc[1]).T.reset_index(drop = True)
    sadf1 = sadfone.add_suffix("_1").drop(['SAMPLE_ALIAS_1', 'LIBRARY_1', 'WORST_CXT_1', 
                                       'WORST_CXT_QSCORE_1', 'WORST_PRE_CXT_1', 
                                       'WORST_PRE_CXT_QSCORE_1', 'WORST_POST_CXT_1',
                                       'WORST_POST_CXT_QSCORE_1'], axis = 1)
    sadftwo = pd.DataFrame(table.loc[2]).T.reset_index(drop = True)
    sadf2 = sadftwo.add_suffix("_2").drop(['SAMPLE_ALIAS_2', 'LIBRARY_2', 'WORST_CXT_2', 
                                       'WORST_CXT_QSCORE_2', 'WORST_PRE_CXT_2', 
                                       'WORST_PRE_CXT_QSCORE_2', 'WORST_POST_CXT_2',
                                       'WORST_POST_CXT_QSCORE_2'], axis = 1)
    sadfthree = pd.DataFrame(table.loc[3]).T.reset_index(drop = True)
    sadf3 = sadfthree.add_suffix("_3").drop(['SAMPLE_ALIAS_3', 'LIBRARY_3','WORST_CXT_3', 
                                        'WORST_CXT_QSCORE_3','WORST_PRE_CXT_3', 
                                         'WORST_PRE_CXT_QSCORE_3', 'WORST_POST_CXT_3',
                                       'WORST_POST_CXT_QSCORE_3'], axis = 1)
    sadffour = pd.DataFrame(table.loc[4]).T.reset_index(drop = True)
    sadf4 = sadffour.add_suffix("_4").drop(['SAMPLE_ALIAS_4', 'LIBRARY_4', 'WORST_CXT_4', 
                                        'WORST_CXT_QSCORE_4', 'WORST_PRE_CXT_4', 
                                        'WORST_PRE_CXT_QSCORE_4', 'WORST_POST_CXT_4',
                                       'WORST_POST_CXT_QSCORE_4'], axis = 1)
    sadffive = pd.DataFrame(table.loc[5]).T.reset_index(drop = True)
    sadf5 = sadffive.add_suffix("_5").drop(['SAMPLE_ALIAS_5', 'LIBRARY_5', 
                                        'WORST_CXT_5', 'WORST_CXT_QSCORE_5', 'WORST_PRE_CXT_5', 
                                        'WORST_PRE_CXT_QSCORE_5', 'WORST_POST_CXT_5',
                                       'WORST_POST_CXT_QSCORE_5'], axis = 1)
    sadfsix = pd.DataFrame(table.loc[6]).T.reset_index(drop = True)
    sadf6 = sadfsix.add_suffix("_6").drop(['SAMPLE_ALIAS_6', 'LIBRARY_6', 
                                       'WORST_CXT_6', 'WORST_CXT_QSCORE_6', 
                                       'WORST_PRE_CXT_6', 'WORST_PRE_CXT_QSCORE_6', 'WORST_POST_CXT_6',
                                       'WORST_POST_CXT_QSCORE_6'], axis = 1)
    sadfseven = pd.DataFrame(table.loc[7]).T.reset_index(drop = True)
    sadf7 = sadfseven.add_suffix("_7").drop(['SAMPLE_ALIAS_7', 'LIBRARY_7', 
                                         'WORST_CXT_7', 'WORST_CXT_QSCORE_7', 
                                         'WORST_PRE_CXT_7', 'WORST_PRE_CXT_QSCORE_7', 'WORST_POST_CXT_7',
                                       'WORST_POST_CXT_QSCORE_7'], axis = 1)
    sadfeight = pd.DataFrame(table.loc[8]).T.reset_index(drop = True)
    sadf8 = sadfeight.add_suffix("_8").drop(['SAMPLE_ALIAS_8', 'LIBRARY_8', 'WORST_CXT_8', 
                                         'WORST_CXT_QSCORE_8', 'WORST_PRE_CXT_8', 
                                         'WORST_PRE_CXT_QSCORE_8', 'WORST_POST_CXT_8',
                                       'WORST_POST_CXT_QSCORE_8'], axis = 1)
    sadfnine = pd.DataFrame(table.loc[9]).T.reset_index(drop = True)
    sadf9 = sadfnine.add_suffix("_9").drop(['SAMPLE_ALIAS_9', 'LIBRARY_9', 
                                        'WORST_CXT_9', 'WORST_CXT_QSCORE_9', 
                                        'WORST_PRE_CXT_9', 'WORST_PRE_CXT_QSCORE_9', 'WORST_POST_CXT_9',
                                       'WORST_POST_CXT_QSCORE_9'], axis = 1)
    sadften = pd.DataFrame(table.loc[10]).T.reset_index(drop = True)
    sadf10 = sadften.add_suffix("_10").drop(['SAMPLE_ALIAS_10', 'LIBRARY_10', 
                                         'WORST_CXT_10', 'WORST_CXT_QSCORE_10', 
                                         'WORST_PRE_CXT_10', 'WORST_PRE_CXT_QSCORE_10', 'WORST_POST_CXT_10',
                                       'WORST_POST_CXT_QSCORE_10'], axis = 1)
    sadfeleven = pd.DataFrame(table.loc[11]).T.reset_index(drop = True)
    sadf11 = sadfeleven.add_suffix("_11").drop(['SAMPLE_ALIAS_11', 'LIBRARY_11', 
                                            'WORST_CXT_11', 'WORST_CXT_QSCORE_11', 
                                            'WORST_PRE_CXT_11', 'WORST_PRE_CXT_QSCORE_11', 'WORST_POST_CXT_11',
                                       'WORST_POST_CXT_QSCORE_11'], axis = 1)
    sadftwelve = pd.DataFrame(table.loc[12]).T.reset_index(drop = True)
    sadf12 = sadftwelve.add_suffix("_12").drop(['SAMPLE_ALIAS_12', 'LIBRARY_12', 
                                            'WORST_CXT_12', 'WORST_CXT_QSCORE_12', 'WORST_PRE_CXT_12', 
                                            'WORST_PRE_CXT_QSCORE_12', 'WORST_POST_CXT_12',
                                       'WORST_POST_CXT_QSCORE_12'], axis = 1)
    SAMresult = pd.concat([sadf1,sadf2, sadf3, sadf4, sadf5, sadf6, sadf7, sadf8, sadf9, sadf10, sadf11, sadf12 ], axis = 1)
    return SAMresult

def windows(table, window):
    if len(table.columns) == 2:
        return window_metrics_mean(table, window)
    
    if len(table.columns) == 7:
        return window_metrics_std(table, window)
        

def mainfile(index):
    allsamplecolumns = ['alignment_summary_metrics', 'base_distribution_by_cycle_table', 'gc_bias_summary_metrics', 
                        'insert_size_metrics', 'mean_quality_by_cycle_table', 
                        'sequencing_artifact_summary_metrics', 'quality_score_table', 'quality_yield_metrics', 
                        'raw_wgs_metrics', 'wgs_metrics']
    table_name = "sample"
    samples = pd.read_csv(io.StringIO(fiss.fapi.get_entities_tsv(project, workspace, 'sample').text), sep='\t')
    samples.rename(columns = {'entity:sample_id':'sample'}, inplace = True)
#     specificcolumns = samples[['alignment_summary_metrics', 'base_distribution_by_cycle_table', 'gc_bias_summary_metrics', 
#                         'insert_size_metrics', 'mean_quality_by_cycle_table', 
#                         'sequencing_artifact_summary_metrics', 'quality_score_table', 'quality_yield_metrics', 
#                         'raw_wgs_metrics', 'wgs_metrics']]
#     dropemptyrows = specificcolumns.dropna(axis = 0)
#     dropemptycolumns = dropemptyrows.dropna(axis = 1)
#     files = ! ls
    
    Dict = {}
    for i in index:
        result = []
        for column in allsamplecolumns:
            metric = samples.at[i, column]
            if type(metric) == float and math.isnan(metric):
                continue
            else:
                result.append(metric)
        if result != []:
            key1 = result[0].split('.')[0].split('/')[-1]
            Dict[key1] = result
    return Dict
            
#     for column in allsamplecolumns:
#         for i in index:
#             metric = dropemptycolumns[column][i]
#             if metric == 'NaN':
#                 continue
#             else:
#                 ! gsutil cp $metric .
#                 result.append(metric.split("/")[-1])
#         return result
        
    
sampledict = mainfile(range(50))
allsamplecolumns = ['alignment_summary_metrics', 'base_distribution_by_cycle_table', 'gc_bias_summary_metrics', 
                        'insert_size_metrics', 'mean_quality_by_cycle_table', 
                        'sequencing_artifact_summary_metrics', 'quality_score_table', 'quality_yield_metrics', 
                        'raw_wgs_metrics', 'wgs_metrics']
def concattables():
    result = pd.DataFrame()
    for k in sampledict.keys():
        listoftables = []
        for v in sampledict.get(k):
            if ('alignment_summary_metrics' in v): 
                alignment1 = addrows(v, '## METRICS CLASS')
                suffixes = alignment1['CATEGORY'].tolist()
                alignment = concatrows(alignment1, suffixes)
                listoftables.append(alignment)
            elif ('base_distribution_by_cycle_table' in v):
                base = windows(addrows(v, '## METRICS CLASS'), 50)
                listoftables.append(base)
            elif('gc_bias_summary_metrics' in v):
                gc = addrows(v, '## METRICS CLASS')
                gcdrop = gc.reset_index(drop = True)
                listoftables.append(gcdrop)
            elif ('insert_size_metrics' in v):
                insert1 = addrows(v, '## METRICS CLASS')
                suffixes = insert1['PAIR_ORIENTATION'].tolist()
                insert = concatrows(insert1, suffixes)
                listoftables.append(insert)
            elif('mean_quality_by_cycle_table' in v):
                meanquality = windows(addrows(v, '## HISTOGRAM'), 50)
                listoftables.append(meanquality)
            elif('pre_adapter_summary_metrics' in v):
                sequencing = sequencing_artifact_metrics(addrows(v, '## METRICS CLASS'))
                listoftables.append(sequencing)
            elif('quality_score_distribution' in v):
                qs1 = addrows(v, '## HISTOGRAM')
                suffixes = qs1['QUALITY'].tolist()
                qs = concatrows(qs1, suffixes)
                listoftables.append(qs)
            elif('quality_yield_metrics' in v):
                qym = addrows(v, '## METRICS CLASS')
                qymdrop = qym.reset_index(drop = True)
                listoftables.append(qymdrop)
            elif('raw_wgs_metrics' in v):
                raw = addrows(v, '## METRICS CLASS')
                rawdrop = raw.reset_index(drop = True).add_suffix('raw')
                listoftables.append(rawdrop)
            elif('wgs_metrics' in v):
                wgs = addrows(v, '## METRICS CLASS')
                wgsdrop = wgs.reset_index(drop = True)
                listoftables.append(wgsdrop)
        newrow = pd.concat(listoftables, axis = 1)
        test =[]
#         for i in range(len(newrow.columns)):
#             for j in range(i + 1, len(newrow.columns)):
#                 if newrow.columns[i] == newrow.columns[j]:
#                     test.append(newrow.columns[i])
#         print(len(test))
#         if len(set(newrow.columns)) != len(newrow.columns):
#             print('duplicate')
        
                
#         duplicates = find_duplicates(newrow)
#         print(duplicates)
        result = result.loc[:, ~result.columns.duplicated()].copy()
        newrow = newrow.loc[:, ~newrow.columns.duplicated()].copy()
        result = pd.concat([result, newrow], ignore_index = True)
    return result
        

        