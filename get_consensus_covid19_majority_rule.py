#! /powerapps/share/python-anaconda-3.6/bin/python

from optparse import OptionParser
import glob
import tqdm
import pandas as pd
from consensus_tools import estimate_insertion_freq, post_processing, deletion_group_filter, insertion_group_filter


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with freqs files")
    parser.add_option('-o', '--o', dest='output_basename')
    parser.add_option('-m', '--mask', action='store_true', default=False, dest='mask', help='whether to mask streches of known bases shorter than 20 amongst stretches of Ns')

    (options, args) = parser.parse_args()
    directory = options.directory
    output_basename = options.output_basename
    freqs_files = glob.glob(f"{directory}/*freqs.csv")
    con_all = ""
    mutations = pd.DataFrame()
    mutations_turned_ns = pd.DataFrame()
    missing_poss = pd.DataFrame()
    end_pos = 29903

    for f in tqdm.tqdm(freqs_files):
        if "all.freqs.csv" in f or 'all.csv' in f or 'all_freqs.csv' in f:
            continue
        sample = f.split("/")[-1].split("_")[0]

        df = pd.read_csv(f, '\t')
        df = df[(df.Read_count > 10)]
        df_snps = df[
            (df["Rank"] == 0) & (df["Ref"] != "-") & (df["Pos"] >= 55) & (df["Pos"] <= 29836) & (
                        df["Base"] != '-')]
        df_deletions = deletion_group_filter(df, initial_threshold=0.75 * 0.5, final_threshold=0.5)
        df = estimate_insertion_freq(df)
        df_insertions = insertion_group_filter(df, initial_threshold=0.5 * 0.75, final_threshold=0.5)
        df_changes = pd.concat([df_snps, df_deletions, df_insertions], sort=False).sort_values("Pos")

        count = int(df_changes.head(1)["Pos"])
        loc = int(df_changes.head(1)["Pos"]) - 1
        con = (count) * "N"
        for index, row in df_changes.iterrows():
            new_loc = int(row["Pos"])
            coverage = int(row["Read_count"])
            ref_base = str(row["Ref"])
            base = str(row["Base"])
            frequency = float(row["Freq"])

            if ref_base != '-':  # not insertion
                if new_loc - loc != 1:
                    con += "N" * (new_loc - loc - 1)
                    missing_poss = missing_poss.append({"start": loc, "end": new_loc, "sample": sample},
                                                       ignore_index=True)

                if base == ref_base:
                    con += base
                else:  # Rank0 ref_base != base
                    con += base
                    mutations = mutations.append(
                            {"Pos": new_loc, "Ref": ref_base, "Base": base, "Read_count": coverage,
                             "Freq": frequency, "sample": sample}, ignore_index=True)
                loc = new_loc

            if ref_base == '-' and new_loc - loc < 1:  # insertion
                con += base
                mutations = mutations.append(
                    {"Pos": row["Pos"], "Ref": ref_base, "Base": base, "Read_count": coverage,
                     "Freq": frequency, "sample": sample}, ignore_index=True)

        con += (end_pos - loc - 1) * "N"
        if options.mask:
            con = post_processing(con)
        con_all += f">Israel/{sample}/2021\n{con}\n"

    with open(f"{directory}/{output_basename}_consensus_all.fasta", "w") as handle:
        handle.write(con_all)
    missing_poss.to_csv(f"{directory}/{output_basename}_missing_positions_all.csv", index=False)
    mutations.to_csv(f"{directory}/{output_basename}_mutations_all.csv", index=False)
    mutations_turned_ns.to_csv(f"{directory}/{output_basename}_mutations_turned_ns_all.csv", index=False)

if __name__ == "__main__":
    main()

