#! /powerapps/share/python-anaconda-3.6/bin/python

import pandas as pd
import re

def estimate_insertion_freq(df, extra_columns=[]):
    '''
    This function gets a freqs file(s) dataframe, calculates the frequency of insertions by using the read count of the
    previous base and returns a dataframe including this.
    :param df: a dataframe of freqs file(s).
    :param extra_columns: if df contains more than the basic freqs columns, for example a column of Sample_id etc.,
    provide a list of the extra columns to be included.
    :return: df with extra columns describing insertion frequency.
    '''
    read_counts = df[(df.Ref != '-')][ extra_columns + ["Pos", "Read_count"]].drop_duplicates()
    read_counts.rename(columns={"Read_count":'estimated_read_count', "Pos":'rounded_pos'}, inplace=True)
    insertions = df[(df.Ref == '-')]
    not_insertions = df[(df.Ref != '-')]
    insertions['rounded_pos'] = insertions.Pos.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on= extra_columns + ['rounded_pos'])
    insertions['estimated_freq'] = insertions.Freq * insertions.Read_count / insertions.estimated_read_count
    df = pd.concat([insertions, not_insertions], sort=False)
    return df.sort_values(extra_columns + ["Pos"])

def post_processing(con_seq):
    def repl(m):
        return 'N' * len(m.group())

    return re.sub('N[AGCT-]{1,20}N', repl, (re.sub('N[AGCT-]{1,20}N', repl, con_seq)))

def ranges(nums):
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))

def deletion_group_filter(df, initial_threshold, final_threshold):
    deletions_df = df[(df["Freq"] >= initial_threshold) & (df["Ref"] != "-") & (df["Pos"] >= 55) & (df["Pos"] <= 29836) & (df["Base"] == '-')]
    deletions_df_to_keep = []
    rs = ranges(deletions_df.sort_values("Pos").Pos.tolist())
    for r in rs:
        r_df = deletions_df[(deletions_df.Pos >= r[0]) & (deletions_df.Pos <= r[1])]
        if r_df.Freq.mean() >= final_threshold:
            deletions_df_to_keep.append(r_df)
    if deletions_df_to_keep != []:
        return pd.concat(deletions_df_to_keep)
    return pd.DataFrame()

def insertion_group_filter(df, initial_threshold, final_threshold):
    insertions_df = df[(df["estimated_freq"] >= initial_threshold) & (df["Ref"] == "-") & (df["Base"] != "-") & (df["Pos"] >= 55) & (df["Pos"] <= 29836) & (df['estimated_read_count'] > 10)]
    insertions_df_to_keep = []
    for r in insertions_df.rounded_pos.drop_duplicates().sort_values().tolist():
        r_df = insertions_df[(insertions_df.rounded_pos == r)]
        if r_df.estimated_freq.mean() >= final_threshold:
            insertions_df_to_keep.append(r_df)
    if insertions_df_to_keep != []:
        return pd.concat(insertions_df_to_keep)
    return pd.DataFrame()