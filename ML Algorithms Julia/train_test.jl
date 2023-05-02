"A given dataframe is split into train:test"
function train_test_split(df,df_cont, labels, split_nr; total_split = 500)
    N_raw, n = size(df) #rows and nr predictors of the given data
    N_test, n_cont = size(df_cont)
    N_raw = max(N_raw, N_test) #in case one is 0
    if split_nr > total_split #is not possible
        error("split number is larger than the total allowed split!")
    end
    Random.seed!(split_nr) #set the seed according to the split number
    test_instances = Random.shuffle(1:N_raw)[1:floor(Int, 0.2 * N_raw)]

    #training/test sets
    labels_train = labels[Not(test_instances),:]
    labels_test = labels[test_instances,:]
    if n > 0
        df_train = df[Not(test_instances), :]
        df_test = df[test_instances, :]
    else
        df_train = DataFrame([])
        df_test = DataFrame([])
    end

    if n_cont > 0
        df_cont_train = df_cont[Not(test_instances), :]
        df_cont_test = df_cont[test_instances, :]
    else
        df_cont_train = DataFrame([])
        df_cont_test = DataFrame([])
    end

    return df_train,df_cont_train, labels_train, df_test, df_cont_test, labels_test
end
"A given dataframe is split into train:test"
function train_test_split_simple(df, labels, split_nr; total_split = 500)
    if split_nr > total_split #is not possible
        error("split number is larger than the total allowed split!")
    end
    Random.seed!(split_nr) #set the seed according to the split number
    test_instances = Random.shuffle(1:size(df)[1])[1:floor(Int, 0.2 * size(df)[1])]

    #training/test sets
    labels_train = labels[Not(test_instances),:]
    labels_test = labels[test_instances,:]

    df_train = df[Not(test_instances), :]
    df_test = df[test_instances, :]

    return df_train,df_test, labels_train, labels_test
end