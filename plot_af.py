import matplotlib.pyplot as plt
import argparse
import scipy.stats as sps
import numpy as np
from scipy.interpolate import LSQUnivariateSpline


''' Example Usage: python3 final_plot_af.py -csv Chr_4_allele_frequencies.csv -contig 4 -window_size 100000 -o Chr4.png '''


def parse_csv_line(csv_line):
    contents = csv_line.split(',')
    position = int(contents[1])
    v48_count = int(contents[3])
    v47_count = int(contents[5])
    return position, v48_count, v47_count


def automate_bins(contig_len, window_size):
    num_bins = int(np.ceil((contig_len / window_size)))
    marg_count = np.zeros(shape=num_bins, dtype=np.float128)
    total_count = np.zeros(shape=num_bins, dtype=np.float128)
    bins = np.zeros(shape=num_bins, dtype=np.int32)
    for i in range(num_bins):
        bins[i] = i * window_size
    return marg_count, total_count, bins, num_bins


def plot_af(cargs):
    expected_dict = {'2': 0.25, '3': 0.25, '4': 0.25, '5': 0.25, '6': 0.25, 'X': 1 / 3}
    contig_lens_dict = {'2': 37e6, '3': 29e6, '4': 30e6, '5': 28e6, '6': 22e6, 'X': 32e6}
    contig = cargs.contig
    contig_len = contig_lens_dict[contig]
    expectation = expected_dict[contig]
    window_size = cargs.window_size

    marg_count, total_count, bins, num_bins = automate_bins(contig_len, window_size)
    csv_file = cargs.csv
    output_file = cargs.out_file
    sites = np.zeros(shape=num_bins, dtype=np.int32)
    use_for_read_depth = np.full(shape=num_bins, dtype=bool, fill_value=False)
    fwer = cargs.fwer
    BonferroniP = fwer / num_bins
    threshold = 1 - BonferroniP
    BonferroniCV = sps.norm.ppf(q=threshold)

    with open(csv_file, "r") as csvf:
        next(csvf)
        for csv_line in csvf:
            position, v48_count, v47_count = parse_csv_line(csv_line)
            index = position // window_size
            marg_count[index] += v48_count
            total_count[index] += (v47_count + v48_count)
            sites[index] += 1

    positions = []
    ratios = []
    colors = []
    for i in range(num_bins):
        m_i = marg_count[i]
        t_i = total_count[i]
        n_i = sites[i]
        if (t_i != 0) and (m_i != 0) and (m_i != t_i):
            use_for_read_depth[i] = 1
            ratio = m_i / t_i
            ratios.append(ratio)
            positions.append(bins[i])
            sd_i = np.sqrt((expectation * (1 - expectation)) / n_i)
            Z = np.abs((ratio - expectation) / sd_i)
            if Z > BonferroniCV:
                colors.append(1)
            else:
                colors.append(0)

    spline_flag = cargs.spline
    if spline_flag:
        spline = LSQUnivariateSpline(positions, ratios, k=5, t=positions[1:-1][0::16])
        plt.figure(figsize=(15, 12))
        plt.scatter(positions, ratios, c=colors, cmap='copper', edgecolors='k', s=40)
        plt.plot(positions, spline(positions), color='b', linewidth=4, linestyle='dashed',
                 label='Quintic Spline')
        plt.axhline(y=expectation, color='r', linestyle='-', linewidth=2, label='Expectation')
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.ylabel("Frequency of Vir48 Allele")
        plt.xlabel(f"Positions on Chromosome {contig} ({window_size} bp bins)")
        plt.title(f"Chromosome {contig}'s Plot")
        plt.savefig("Chr" + "_" + str(contig) + "_" + output_file, dpi=100)
    else:
        plt.figure(figsize=(15, 12))
        plt.scatter(positions, ratios, c=colors, cmap='copper', edgecolors='k', s=40)
        plt.axhline(y=expectation, color='r', linestyle='-', linewidth=2, label='Expectation')
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.ylabel("Frequency of Vir48 Allele")
        plt.xlabel(f"Positions on Chromosome {contig} ({window_size} bp bins)")
        plt.title(f"Chromosome {contig}'s Plot")
        plt.savefig("Chr" + "_" + str(contig) + "_" + output_file, dpi=100)

    read_depth_flag = cargs.read_depth
    if read_depth_flag:
        rd_avgs = []
        for i in range(num_bins):
            if use_for_read_depth[i]:
                avg_i = total_count[i] / sites[i]
                rd_avgs.append(avg_i)
        plt.figure(figsize=(15, 12))
        plt.plot(positions, rd_avgs, color='b')
        plt.ylabel("Average Read Depth")
        plt.xlabel(f"Positions on Chromosome {contig} ({window_size} bp bins)")
        plt.title(f"Chromosome {contig}'s Read Depth Plot")
        plt.savefig("Chr" + "_" + str(contig) + "_" + "read_depth.png", dpi=100)


def main():
    parser = argparse.ArgumentParser(
        description="Plot allele frequencies at bi-allelic sites in Vir47 and Vir48 genomes."
                    "\nWill output plots of the average read depth of bins.")
    parser.add_argument("-csv", action="store", dest="csv", type=str)
    parser.add_argument("-contig", action="store", dest="contig", type=str)
    parser.add_argument("-window_size", action="store", dest="window_size", type=int, default=125000)
    parser.add_argument("-fwer", action="store", dest="fwer", type=float, default=0.01)
    parser.add_argument("-read_depth", action="store", dest="read_depth", type=bool, default=False)
    parser.add_argument("-spline", action="store", dest="spline", type=bool, default=False)
    parser.add_argument("-o", action="store", dest="out_file", type=str,
                        default="allele_frequencies.png")
    args = parser.parse_args()
    plot_af(args)
    return


if __name__ == "__main__":
    main()