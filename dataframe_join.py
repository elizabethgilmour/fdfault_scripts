import os
RESULT = os.path.join("data", "patch_lengthmagdec3.txt")
PROBLEM = os.path.join("problems", "patch_length.txt")

import pandas as pd

df_results = pd.read_csv(RESULT)

df_problems =pd.read_csv(PROBLEM)

df_problems['seed_file'] = list(df_problems.index)

joined = pd.merge(df_results, df_problems, left_on = 'seed', right_on = 'seed_file')

del joined['seed_x']

joined.to_csv(os.path.join('data', 'clean_data', 'sun_dec_3.txt'))