def read_in(file_text, label):
    started = False
    newlist = []
    
    for line in file_text:
        if line.startswith(label):
            started = True
            continue
        if started:
            if line != '\n':
                newlist.append(line.rstrip('\n'))
        if started and line == '\n':
            return newlist
            break
                
def newtable(file_text, label):
    temp = read_in(file_text, label)
    columnnames = temp[0].split("\t")
    values = temp[1].split("\t")
    #values = [float(n) for n in values if n.isnumeric()]
    dataframe = pd.DataFrame(columnnames).T
    dataframe.columns = columnnames
    dataframe_length = len(dataframe)
    dataframe.loc[dataframe_length] = values
    dataframe = dataframe.drop(0)
    return dataframe

def addrows(file_name, label):
    with open(file_name) as file:
        file_text = [line for line in file]
        temp = read_in(file_text, label)
        onerow = newtable(file_text, label)
    for i in range(len(temp)):
        nextrow = temp[i].split("\t")
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
    sadf1 = sadfone.add_suffix("_AC").drop(['SAMPLE_ALIAS_AC', 'LIBRARY_AC', 'WORST_CXT_AC', 
                                       'WORST_CXT_QSCORE_AC', 'WORST_PRE_CXT_AC', 
                                       'WORST_PRE_CXT_QSCORE_AC', 'WORST_POST_CXT_AC',
                                       'WORST_POST_CXT_QSCORE_AC'], axis = 1)
    sadftwo = pd.DataFrame(table.loc[2]).T.reset_index(drop = True)
    sadf2 = sadftwo.add_suffix("_AG").drop(['SAMPLE_ALIAS_AG', 'LIBRARY_AG', 'WORST_CXT_AG', 
                                       'WORST_CXT_QSCORE_AG', 'WORST_PRE_CXT_AG', 
                                       'WORST_PRE_CXT_QSCORE_AG', 'WORST_POST_CXT_AG',
                                       'WORST_POST_CXT_QSCORE_AG'], axis = 1)
    sadfthree = pd.DataFrame(table.loc[3]).T.reset_index(drop = True)
    sadf3 = sadfthree.add_suffix("_AT").drop(['SAMPLE_ALIAS_AT', 'LIBRARY_AT','WORST_CXT_AT', 
                                        'WORST_CXT_QSCORE_AT','WORST_PRE_CXT_AT', 
                                         'WORST_PRE_CXT_QSCORE_AT', 'WORST_POST_CXT_AT',
                                       'WORST_POST_CXT_QSCORE_AT'], axis = 1)
    sadffour = pd.DataFrame(table.loc[4]).T.reset_index(drop = True)
    sadf4 = sadffour.add_suffix("_CA").drop(['SAMPLE_ALIAS_CA', 'LIBRARY_CA', 'WORST_CXT_CA', 
                                        'WORST_CXT_QSCORE_CA', 'WORST_PRE_CXT_CA', 
                                        'WORST_PRE_CXT_QSCORE_CA', 'WORST_POST_CXT_CA',
                                       'WORST_POST_CXT_QSCORE_CA'], axis = 1)
    sadffive = pd.DataFrame(table.loc[5]).T.reset_index(drop = True)
    sadf5 = sadffive.add_suffix("_CG").drop(['SAMPLE_ALIAS_CG', 'LIBRARY_CG', 
                                        'WORST_CXT_CG', 'WORST_CXT_QSCORE_CG', 'WORST_PRE_CXT_CG', 
                                        'WORST_PRE_CXT_QSCORE_CG', 'WORST_POST_CXT_CG',
                                       'WORST_POST_CXT_QSCORE_CG'], axis = 1)
    sadfsix = pd.DataFrame(table.loc[6]).T.reset_index(drop = True)
    sadf6 = sadfsix.add_suffix("_CT").drop(['SAMPLE_ALIAS_CT', 'LIBRARY_CT', 
                                       'WORST_CXT_CT', 'WORST_CXT_QSCORE_CT', 
                                       'WORST_PRE_CXT_CT', 'WORST_PRE_CXT_QSCORE_CT', 'WORST_POST_CXT_CT',
                                       'WORST_POST_CXT_QSCORE_CT'], axis = 1)
    sadfseven = pd.DataFrame(table.loc[7]).T.reset_index(drop = True)
    sadf7 = sadfseven.add_suffix("_GA").drop(['SAMPLE_ALIAS_GA', 'LIBRARY_GA', 
                                         'WORST_CXT_GA', 'WORST_CXT_QSCORE_GA', 
                                         'WORST_PRE_CXT_GA', 'WORST_PRE_CXT_QSCORE_GA', 'WORST_POST_CXT_GA',
                                       'WORST_POST_CXT_QSCORE_GA'], axis = 1)
    sadfeight = pd.DataFrame(table.loc[8]).T.reset_index(drop = True)
    sadf8 = sadfeight.add_suffix("_GC").drop(['SAMPLE_ALIAS_GC', 'LIBRARY_GC', 'WORST_CXT_GC', 
                                         'WORST_CXT_QSCORE_GC', 'WORST_PRE_CXT_GC', 
                                         'WORST_PRE_CXT_QSCORE_GC', 'WORST_POST_CXT_GC',
                                       'WORST_POST_CXT_QSCORE_GC'], axis = 1)
    sadfnine = pd.DataFrame(table.loc[9]).T.reset_index(drop = True)
    sadf9 = sadfnine.add_suffix("_GT").drop(['SAMPLE_ALIAS_GT', 'LIBRARY_GT', 
                                        'WORST_CXT_GT', 'WORST_CXT_QSCORE_GT', 
                                        'WORST_PRE_CXT_GT', 'WORST_PRE_CXT_QSCORE_GT', 'WORST_POST_CXT_GT',
                                       'WORST_POST_CXT_QSCORE_GT'], axis = 1)
    sadften = pd.DataFrame(table.loc[10]).T.reset_index(drop = True)
    sadf10 = sadften.add_suffix("_TA").drop(['SAMPLE_ALIAS_TA', 'LIBRARY_TA', 
                                         'WORST_CXT_TA', 'WORST_CXT_QSCORE_TA', 
                                         'WORST_PRE_CXT_TA', 'WORST_PRE_CXT_QSCORE_TA', 'WORST_POST_CXT_TA',
                                       'WORST_POST_CXT_QSCORE_TA'], axis = 1)
    sadfeleven = pd.DataFrame(table.loc[11]).T.reset_index(drop = True)
    sadf11 = sadfeleven.add_suffix("_TC").drop(['SAMPLE_ALIAS_TC', 'LIBRARY_TC', 
                                            'WORST_CXT_TC', 'WORST_CXT_QSCORE_TC', 
                                            'WORST_PRE_CXT_TC', 'WORST_PRE_CXT_QSCORE_TC', 'WORST_POST_CXT_TC',
                                       'WORST_POST_CXT_QSCORE_TC'], axis = 1)
    sadftwelve = pd.DataFrame(table.loc[12]).T.reset_index(drop = True)
    sadf12 = sadftwelve.add_suffix("_TG").drop(['SAMPLE_ALIAS_TG', 'LIBRARY_TG', 
                                            'WORST_CXT_TG', 'WORST_CXT_QSCORE_TG', 'WORST_PRE_CXT_TG', 
                                            'WORST_PRE_CXT_QSCORE_TG', 'WORST_POST_CXT_TG',
                                       'WORST_POST_CXT_QSCORE_TG'], axis = 1)
    SAMresult = pd.concat([sadf1,sadf2, sadf3, sadf4, sadf5, sadf6, sadf7, sadf8, sadf9, sadf10, sadf11, sadf12 ], axis = 1)
    return SAMresult

def windows(table, window):
    if len(table.columns) == 2:
        return window_metrics_mean(table, window)
    
    if len(table.columns) == 7:
        return window_metrics_std(table, window)
        

def mainfile(index, localize = True):
    allsamplecolumns = ['alignment_summary_metrics', 'base_distribution_by_cycle_table', 'gc_bias_summary_metrics', 
                        'insert_size_metrics', 'mean_quality_by_cycle_table', 
                        'sequencing_artifact_summary_metrics', 'quality_score_table', 'quality_yield_metrics', 
                        'raw_wgs_metrics', 'wgs_metrics']
    table_name = "sample"
    samples = pd.read_csv(io.StringIO(fiss.fapi.get_entities_tsv(project, workspace, 'sample').text), sep='\t')
    samples.rename(columns = {'entity:sample_id':'sample'}, inplace = True)

    ! rm -f filelist.txt
    Dict = {}
    for i in index:
        result = []
        for column in allsamplecolumns:
            metric = samples.at[i, column]
            if type(metric) == float and math.isnan(metric):
                continue
            else:
                if localize and not os.path.exists(metric.split('/')[-1]):
                    ! echo $metric >> filelist.txt
                result.append(metric)
        if result != []:
            key1 = result[0].split('.')[0].split('/')[-1]
            Dict[key1] = result
    return Dict
    

def concattables():
    result = pd.DataFrame()
    for k in sampledict.keys():
        listoftables = []
        for v in sampledict.get(k):
            if ('alignment_summary_metrics' in v):
                alignment1 = addrows((v.split('/')[-1]), '## METRICS CLASS')
                suffixes = alignment1['CATEGORY'].tolist()
                alignment = concatrows(alignment1, suffixes)
                listoftables.append(alignment)
            elif ('base_distribution_by_cycle_table' in v):
                base = windows(addrows((v.split('/')[-1]), '## METRICS CLASS'), 50)
                listoftables.append(base)
            elif('gc_bias_summary_metrics' in v):
                gc = addrows((v.split('/')[-1]), '## METRICS CLASS')
                gcdrop = gc.reset_index(drop = True)
                listoftables.append(gcdrop)
            elif ('insert_size_metrics' in v):
                insert1 = addrows((v.split('/')[-1]), '## METRICS CLASS')
                suffixes = insert1['PAIR_ORIENTATION'].tolist()
                insert = concatrows(insert1, suffixes)
                listoftables.append(insert)
            elif('mean_quality_by_cycle_table' in v):
                meanquality = windows(addrows((v.split('/')[-1]), '## HISTOGRAM'), 50)
                listoftables.append(meanquality)
            elif('pre_adapter_summary_metrics' in v):
                sequencing = sequencing_artifact_metrics(addrows((v.split('/')[-1]), '## METRICS CLASS'))
                listoftables.append(sequencing)
            elif('quality_score_distribution' in v):
                qs1 = addrows((v.split('/')[-1]), '## HISTOGRAM')
                suffixes = qs1['QUALITY'].tolist()
                qs = concatrows(qs1, suffixes)
                listoftables.append(qs)
            elif('quality_yield_metrics' in v):
                qym = addrows((v.split('/')[-1]), '## METRICS CLASS')
                qymdrop = qym.reset_index(drop = True)
                listoftables.append(qymdrop)
            elif('raw_wgs_metrics' in v):
                raw = addrows((v.split('/')[-1]), '## METRICS CLASS')
                rawdrop = raw.reset_index(drop = True).add_suffix('raw')
                listoftables.append(rawdrop)
            elif('wgs_metrics' in v):
                wgs = addrows((v.split('/')[-1]), '## METRICS CLASS')
                wgsdrop = wgs.reset_index(drop = True)
                listoftables.append(wgsdrop)    
        newrow = pd.concat(listoftables, axis = 1)
        newrow['Sample_ID'] = k
        result = result.loc[:, ~result.columns.duplicated()].copy()
        newrow = newrow.loc[:, ~newrow.columns.duplicated()].copy()
        result = pd.concat([result, newrow], ignore_index = True)
    return result

def convert_table(table):
    dropped = table.drop(['CATEGORYFIRST_OF_PAIR', 'CATEGORYSECOND_OF_PAIR', 'CATEGORYPAIR', 'PAIR_ORIENTATIONFR', 'PAIR_ORIENTATIONRF', 'PAIR_ORIENTATIONTANDEM',
                          'MEAN_ALIGNED_READ_LENGTHFIRST_OF_PAIR', 'MEAN_ALIGNED_READ_LENGTHSECOND_OF_PAIR',
                          'MEAN_ALIGNED_READ_LENGTHPAIR', 'FOLD_95_BASE_PENALTY', 'Sample_ID', 'ACCUMULATION_LEVEL',
                          'READS_USED', 'SAMPLE', 'LIBRARY', 'SAMPLEFR', 'LIBRARYFR', 'READ_GROUPFR', 'SAMPLERF',
                          'LIBRARYRF', 'READ_GROUPRF', 'SAMPLETANDEM', 'LIBRARYTANDEM', 'READ_GROUPTANDEM',
                          'REF_BASE_AC', 'ALT_BASE_AC', 'ARTIFACT_NAME_AC', 'REF_BASE_AG', 'ALT_BASE_AG',
                          'ARTIFACT_NAME_AG', 'REF_BASE_AT', 'ALT_BASE_AT', 'ARTIFACT_NAME_AT', 'REF_BASE_CA', 
                          'ALT_BASE_CA', 'ARTIFACT_NAME_CA', 'REF_BASE_CG', 'ALT_BASE_CG', 'ARTIFACT_NAME_CG', 
                          'REF_BASE_CT', 'ALT_BASE_CT', 'ARTIFACT_NAME_CT', 'REF_BASE_GA', 'ALT_BASE_GA', 'ARTIFACT_NAME_GA',
                          'REF_BASE_GC', 'ALT_BASE_GC', 'ARTIFACT_NAME_GC', 'REF_BASE_GT', 'ALT_BASE_GT', 
                          'ARTIFACT_NAME_GT', 'REF_BASE_TA', 'ALT_BASE_TA', 'ARTIFACT_NAME_TA', 'REF_BASE_TC', 'ALT_BASE_TC',
                          'ARTIFACT_NAME_TC', 'REF_BASE_TG', 'ALT_BASE_TG', 'ARTIFACT_NAME_TG'], axis = 1) 
    dropped = dropped.dropna(axis = 'columns')
    for column in dropped:
        dropped[column] = pd.to_numeric(dropped[column], errors = 'coerce')
    dropped = dropped.dropna(axis = 'columns')
    return dropped
        

table = convert_table(finaltable)
table['PAIRRATIO_BASES:READS_ALIGNED'] = table['PF_HQ_ALIGNED_BASESPAIR']/(table['PF_READS_ALIGNEDPAIR'])
table['DIFF_TOTALBASES_PF_HQ_ALIGNED_BASESPAIR'] = (table['TOTAL_BASES'] - table['PF_HQ_ALIGNED_BASESPAIR'])
table['PF_BASES_RATIO'] = table['PF_BASES']/(table['TOTAL_BASES'])
table['PF_HQ_ALIGNED_Q20_BASESPAIR_RATIO'] = table['PF_HQ_ALIGNED_Q20_BASESPAIR']/(table['TOTAL_BASES'])
table['PF_HQ_ALIGNED_BASESPAIR_RATIO'] = table['PF_HQ_ALIGNED_BASESPAIR']/(table['TOTAL_BASES'])
table['PF_HQ_ALIGNED_Q20_BASESSECOND_OF_PAIR_RATIO'] = table['PF_HQ_ALIGNED_Q20_BASESSECOND_OF_PAIR']/(table['TOTAL_BASES'])
table['COUNT_OF_Q20_RATIO'] = table['COUNT_OF_Q20']/(table['TOTAL_BASES'])
table['PF_HQ_ALIGNED_BASESSECOND_OF_PAIR_RATIO'] = table['PF_HQ_ALIGNED_BASESSECOND_OF_PAIR']/(table['TOTAL_BASES'])
table['Q20_BASES_RATIO'] = table['Q20_BASES']/(table['TOTAL_BASES'])
table['COUNT_OF_Q10_RATIO'] = table['COUNT_OF_Q10']/(table['TOTAL_BASES'])
table['DIFF_TOTALBASES_PF_HQ_ALIGNED_BASESPAIR_RATIO'] = table['DIFF_TOTALBASES_PF_HQ_ALIGNED_BASESPAIR']/(table['TOTAL_BASES'])
table['PF_HQ_ALIGNED_Q20_BASESFIRST_OF_PAIR_RATIO'] = table['PF_HQ_ALIGNED_Q20_BASESFIRST_OF_PAIR']/(table['TOTAL_BASES'])

def y_list():
    y = []
    tolist = finaltable['Sample_ID'].tolist()
    Dict = is_bad()
    for i in tolist:
        if i in Dict.keys():
            y.append(Dict[i])
    return y


table['qc_pass'] = y_list()
noqcpass = table.drop(['qc_pass'], axis = 1)


X = noqcpass
y = y_list()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
clf = LogisticRegression(random_state=0).fit(X_train, y_train)
clfvalues = list(clf.coef_[0])
columns = list(X.columns)

def merge(list1, list2):
    merged_list = tuple(zip(list1, list2))
    merged_list = sorted(merged_list, key = lambda x: abs(x[0]), reverse = True)
    return merged_list
        


def logistic_reg_predict(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
    clf = LogisticRegression(random_state=0).fit(X_train, y_train)
    return clf.predict(X_test)

def logistic_reg_score(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
    clf = LogisticRegression(random_state=0).fit(X_train, y_train)
    return clf.score(X_test, y_test)

def logistic_reg_decision(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
    clf = LogisticRegression(random_state=0).fit(X_train, y_train)
    return clf.decision_function(X_test)

def logistic_reg_prob(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
    clf = LogisticRegression(random_state=0).fit(X_train, y_train)
    return clf.predict_proba(X_test)

def logistic_reg_log_prob(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
    clf = LogisticRegression(random_state=0).fit(X_train, y_train)
    return clf.predict_log_proba(X_test)

def logistic_reg_model(X, y):
    predict = logistic_reg_predict(X, y)
    score = logistic_reg_score(X, y)
    decision = logistic_reg_decision(X, y)
    probability = logistic_reg_prob(X, y)
    logprob = logistic_reg_log_prob(X, y)
    return predict, score, decision, probability, logprob

cv = KFold(n_splits = 10, random_state = 1, shuffle = True)
model = LogisticRegression()
scores = cross_val_score(model, X, y, scoring = 'accuracy', cv=cv, n_jobs=-1)
print(scores)
print('Accuracy: %.3f (%.3f)' % (mean(scores), stdev(scores)))


topten = table[['PF_BASES_RATIO', 'PF_HQ_ALIGNED_Q20_BASESPAIR_RATIO', 'PF_HQ_ALIGNED_BASESPAIR_RATIO', 'PF_HQ_ALIGNED_Q20_BASESSECOND_OF_PAIR_RATIO',
                      'COUNT_OF_Q20_RATIO', 'DIFF_TOTALBASES_PF_HQ_ALIGNED_BASESPAIR_RATIO', 'PF_HQ_ALIGNED_BASESSECOND_OF_PAIR_RATIO', 'Q20_BASES_RATIO', 'COUNT_OF_Q10_RATIO', 'PF_HQ_ALIGNED_Q20_BASESFIRST_OF_PAIR_RATIO', 'qc_pass']]
palette = {
    '0': 'tab:blue',
    '1': 'tab:red'
}
sns.pairplot(data = topten, vars = topten.columns[:-1], hue = "qc_pass", kind = 'scatter', plot_kws =dict(alpha = 0.2), diag_kind='hist')

hist1 = table[['PF_BASES_RATIO', 'qc_pass']]
sns1 = sns.histplot(data = hist1, x = table['PF_BASES_RATIO'], hue = "qc_pass")
hist2 = table[['PF_HQ_ALIGNED_Q20_BASESPAIR_RATIO', 'qc_pass']]
sns2 = sns.histplot(data = hist2, x = table['PF_HQ_ALIGNED_Q20_BASESPAIR_RATIO'], hue = "qc_pass")
hist3 = table[['PF_HQ_ALIGNED_BASESPAIR_RATIO', 'qc_pass']]
sns3 = sns.histplot(data = hist3, x = table['PF_HQ_ALIGNED_BASESPAIR_RATIO'], hue = "qc_pass")
hist4 = table[['PF_HQ_ALIGNED_Q20_BASESSECOND_OF_PAIR_RATIO', 'qc_pass']]
sns4 = sns.histplot(data = hist4, x = table['PF_HQ_ALIGNED_Q20_BASESSECOND_OF_PAIR_RATIO'], hue = "qc_pass")
hist5 = table[['COUNT_OF_Q20_RATIO', 'qc_pass']]
sns5 = sns.histplot(data = hist5, x = table['COUNT_OF_Q20_RATIO'], hue = "qc_pass")
hist6 = table[['DIFF_TOTALBASES_PF_HQ_ALIGNED_BASESPAIR_RATIO', 'qc_pass']]
sns6 = sns.histplot(data = hist6, x = table['DIFF_TOTALBASES_PF_HQ_ALIGNED_BASESPAIR_RATIO'], hue = "qc_pass")
hist7 = table[['PF_HQ_ALIGNED_BASESSECOND_OF_PAIR_RATIO', 'qc_pass']]
sns7 = sns.histplot(data = hist7, x = table['PF_HQ_ALIGNED_BASESSECOND_OF_PAIR_RATIO'], hue = "qc_pass")
hist8 = table[['Q20_BASES_RATIO', 'qc_pass']]
sns8 = sns.histplot(data = hist8, x = table['Q20_BASES_RATIO'], hue = "qc_pass")
hist9 = table[['COUNT_OF_Q10_RATIO', 'qc_pass']]
sns9 = sns.histplot(data = hist9, x = table['COUNT_OF_Q10_RATIO'], hue = "qc_pass")
hist10 = table[['PF_HQ_ALIGNED_Q20_BASESFIRST_OF_PAIR_RATIO', 'qc_pass']]
sns10 = sns.histplot(data = hist10, x = table['PF_HQ_ALIGNED_Q20_BASESFIRST_OF_PAIR_RATIO'], hue = "qc_pass")


rf = RandomForestRegressor(n_estimators = 100)
rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)
y_pred_list = []
for i in y_pred:
    if (i is not 0) or (i is not 1):
        rounded = round(i)
        y_pred_list.append(rounded)

