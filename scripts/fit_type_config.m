switch fit_type
    case '96well'
        index_list = [1:96];
    case 'one_row'
        index_list = [4:8:92];
    case 'one_column'
        index_list = [65:72];
    case 'one_cross'
        index_list = [4:8:60,65:72,76,84,92];
    case 'gluc_gradient'
        index_list = [1:12];
end