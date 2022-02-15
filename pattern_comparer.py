# -*- coding: utf-8 -*-
"""
Script to compare output of the mass_pattern_finder with a reference file

(C) Mathijs Mabesoone, ETH Zurich
February 2022
"""
from collections import defaultdict
import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import cmocean
from numpy import log10
from math import ceil, floor
import numpy as np
from copy import deepcopy

def sum_peaks(peaks):
    print(sum([len(peaks[rt]) for rt in peaks]))

def obtain_peaks_from_file(file):
    #Extract the peaks from the file
    txt = open(file, 'r').read()
    timepoints = txt.split('Found matching pattern at ')[1:]
    peaks_from_file = dict()
    for timepoint in timepoints:
        masses = timepoint.split('\n\tMass: ')
        time = float(masses.pop(0).replace(' min:',''))
        peaks_from_file[time] = list()
        for mass_dat in masses:
            mass = float(mass_dat.split('\t')[0])
            intensity = float(mass_dat.split('\t')[1].split(': ')[1])
            formulas = [formula.replace('\n', '') for formula in mass_dat.split('\t\t')[1:]]
            peaks_from_file[time].append([mass, intensity, formulas])
    return peaks_from_file

def parse_files(files, args, sample):
    # Pare the input files and obtain the peaks. Cluster them if necessary
    peaks = defaultdict(list)
    for file in files:
        peaks[file] = obtain_peaks_from_file(file)
        if len(peaks[file]) == 0:
            sys.stdout.write(f"Did not find any peaks in {file}.\n")
    files = [file for file in files if len(peaks[file]) > 0 ]

    peaks_in_all_files = peaks[files[0]]

    # If there are more files, make sure the peaks are in all the sample files
    for file_iterator, file in enumerate(files):
        if file_iterator == 0:
            sys.stdout.write(f"Found {len(files)} {''.join(['input' if sample else 'reference'])} files:\n\t" + '\n\t'.join([file for file in files]) + '\n')

        for reference_file in files[file_iterator+1:]:
            sys.stdout.write(f"Finding {''.join(['unique' if sample else 'shared'])} peaks between {os.path.basename(file)} and {os.path.basename(reference_file)}.\n")
            if file_iterator == 0:
                unique_peaks, peaks_in_all_files = compare_peaks(peaks[file], peaks[reference_file], args, sample)
            else:
                unique_peaks, peaks_in_all_files = compare_peaks(peaks_in_all_files, peaks[reference_file], args, sample)
            if not sample:
                for rt in unique_peaks:
                    peaks_in_all_files[rt].extend(unique_peaks[rt])
    return peaks_in_all_files


def compare_input_mass_to_reference_masses(input_mass, reference_masses, args, mode):
    # Small function to see if an input mass is present in a list of reference masses
    result = [None, None]
    for reference_mass in reference_masses:
        if abs(input_mass[0]-reference_mass[0]) < args.mass_tolerance*input_mass[0]:
            if mode == 'compare':
                result[0] = input_mass
            else:
                result[0] = [round((input_mass[0] + reference_mass[0])/2, 4), round((input_mass[1] + reference_mass[1])/2, 1), input_mass[2]]
                break
            return result

    # If no matches are found, the peak is unique
    if input_mass[1] > args.min_intensity:
        result[1] = input_mass
    return result


def compare_peaks(input_peaks, reference_peaks, args, mode):
    # Compare the peaks and return a dictionary of unique and mutual peaks
    # For the mutual peaks, the average of the masses, retention times and intensities is reported

    unique_peaks = defaultdict(list)
    mutual_peaks = defaultdict(list)

    num_reference_peaks = len(reference_peaks)
    rt_ref_peaks = [rt for rt in reference_peaks]

    ref_iterator = 0
    total_input_peaks = len(input_peaks)
    for rt_index, rt in enumerate(input_peaks):
        #if rt > 6.46:
        #     break
        # Check for every elution time...
        if mode == 'compare':
            sys.stdout.write(f"\rProgress: {round(rt_index/total_input_peaks*100, 2)}%")
        # Check if the peaks at all elution times in the sample match to a reference
        if ref_iterator < len(reference_peaks)-1:
            # If the reference iterator for the reference is still much earlier, first catch up
            while rt_ref_peaks[ref_iterator]-rt < -1*args.time_tolerance and ref_iterator < len(rt_ref_peaks)-1:
                ref_iterator += 1

        if rt_ref_peaks[ref_iterator]-rt > args.time_tolerance:
            unique_peaks[rt] = input_peaks[rt]
            continue
        if rt_ref_peaks[-1] < rt:
            unique_peaks[rt] = input_peaks[rt]
            continue

        elif abs(rt_ref_peaks[ref_iterator]-rt) <= args.time_tolerance and ref_iterator < len(rt_ref_peaks):
           input_masses = input_peaks[rt]
           # If the rt range of interest is in the reference, look in the neighboring rts for all masses
           reference_masses = reference_peaks[rt_ref_peaks[ref_iterator]]
           ref_mass_iterator = ref_iterator
           while rt_ref_peaks[ref_mass_iterator]-rt < 1*args.time_tolerance and ref_mass_iterator < len(rt_ref_peaks)-1:
               reference_masses.extend(reference_peaks[rt_ref_peaks[ref_mass_iterator]])
               ref_mass_iterator += 1
           # And then compare each input_mass to the masses in the surrounding reference spectra
           for input_mass in input_masses:
               mutual_peak, unique_peak = compare_input_mass_to_reference_masses(input_mass, reference_masses, args, mode)
               if mutual_peak:
                   mutual_peaks[round(rt, 2)].append(mutual_peak)
               elif unique_peak:
                   unique_peaks[round(rt, 2)].append(unique_peak)

    sys.stdout.write('\n')
    return unique_peaks, mutual_peaks

def filter_peaks_for_identical_masses(input_peaks, args):
    # Combine hits that have masses within the mass tolerance

    # First filter each spectrum for identical peaks
    filtered_peaks = defaultdict(list)
    for rt in input_peaks:
        peaks = sorted(input_peaks[rt], key=lambda peak: peak[0])
        peak_mass, peak_intensity, peak_formula = [0, 0, '']
        start_mass = 0
        peak_iterator = 0
        while peak_iterator < len(peaks)-1:
            mass, intensity, formula = peaks[peak_iterator]
            if mass-start_mass <= 2*args.mass_tolerance*start_mass and intensity > peak_intensity:
                peak_mass = mass
                peak_intensity = intensity
                peak_formula = formula
            if mass-start_mass > 2*args.mass_tolerance*start_mass:
                if peak_iterator > 0:
                    filtered_peaks[rt].append([peak_mass, peak_intensity, peak_formula])
                start_mass = mass
                peak_mass = mass
                peak_intensity = intensity
                peak_formula = formula
            peak_iterator += 1

    rts = [rt for rt in filtered_peaks]

    # Make a copy
    centered_peaks = deepcopy(filtered_peaks)

    for rt_index, rt in enumerate(filtered_peaks):
        peaks = filtered_peaks[rt]
        for mass, intensity, formula in peaks:
            # Check for each peak in the neighborhing retention times
            rt_iterator = 1
            while rt_index+rt_iterator < len(rts):
                next_rt = rts[rt_index+rt_iterator]
                if next_rt-rt > args.time_tolerance:
                    break
                next_peaks = filtered_peaks[next_rt]
                for next_mass, next_intensity, next_formula in next_peaks:
                    if next_mass - mass < -1*args.mass_tolerance*mass:
                        continue
                    elif next_mass - mass > args.mass_tolerance*mass:
                        break
                    elif intensity >= next_intensity:
                        try:
                            centered_peaks[next_rt].remove([next_mass, next_intensity, next_formula])
                        except ValueError:
                            pass
                    elif next_intensity > intensity:
                        try:
                            centered_peaks[rt].remove([mass, intensity, formula])
                        except ValueError:
                            pass
                rt_iterator += 1

    return centered_peaks


def plot_results_in_2D(unique_peaks, plot_filename, time_range, mass_range, args):
    # Plot the analyzed spectra in a single graph
    plot_list = list()
    for rt in unique_peaks:
        spectrum = unique_peaks[rt]
        for peak in spectrum:
            plot_list.append([peak[0], rt, log10(peak[1])])
    plot_list = pd.DataFrame(plot_list, columns=['mass', 'time', 'intensity'])
    plot_list = plot_list.sort_values(by='intensity', ascending='True', ignore_index=True)

    if os.path.isdir(args.input):
        sample_name = os.path.abspath(args.input).split('\\')[-1]
    else:
        sample_name = os.path.basename(args.input)
    if os.path.isdir(args.reference):
        reference_name = os.path.abspath(args.reference).split('\\')[-1]
    else:
        reference_name = os.path.basename(args.reference)

    # Define the plot
    zrange = [power for power in range(floor(min(plot_list['intensity'])), ceil(max(plot_list['intensity'])))]
    plt.figure(figsize=(12, 10))
    plt.title(f"Unique masses\nInput: {sample_name}    Ref: {reference_name}", fontsize=16)
    plot = plt.scatter(plot_list['time'], plot_list['mass'], c=plot_list['intensity'],
                       marker='.',
                       cmap=cmocean.cm.rain,
                       vmin=zrange[0], vmax=zrange[-1])
    plt.xlabel('Time (min)', fontsize=18)
    plt.ylabel('m/z', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((min(time_range), max(time_range)))
    plt.ylim((min(mass_range), max(mass_range)))
    plt.grid(which='both', alpha=0.3)
    cbar = plt.colorbar(plot, shrink=0.5)
    cbar.ax.set_ylabel('log(intensity)', rotation=270, labelpad=20, fontsize=18)
    cbar.set_ticks(zrange)
    cbar.ax.tick_params(labelsize=14)

    # Save the plot
    plot_filename = os.path.splitext(plot_filename)[0]
    if not args.overwrite and any([os.path.isfile(plot_filename+'.svg'), os.path.isfile(plot_filename+'.png')]):
        suffix = 1
        while any([os.path.isfile(plot_filename+'_'+str(suffix)+'.svg'), os.path.isfile(plot_filename+'_'+str(suffix)+'.png')]):
            suffix += 1
        plot_filename = f"{plot_filename}_{suffix}"
    plt.savefig(plot_filename+'.svg', transparent=True, dpi=300)
    plt.savefig(plot_filename+'.png', transparent=True, dpi=300)
    sys.stdout.write(f"Saved plot:\n\t{plot_filename}.png\n\t{plot_filename}.svg\n")

def clean_formula(formula):
    # Remove the parts in the chemical formula that have coefficient 0 and clean
    # the string
    elements = formula.split('_')
    clean_formula = list()
    chemical_groups = ['C2H4','C2H2','CH2','NH3', 'O']
    for element in elements:
        # Check if the last value is 0, not 10/20.. and not a custom element
        if element[-1] == '0' and not ':' in element:
            if element[-2].isdigit() and not any([alkyl in element[-2-len(alkyl):-1] for alkyl in chemical_groups]):
                pass
            else:
                continue
        else:
            clean_formula.append(element)
    customs = [element for element in clean_formula if ':' in element]
    adducts = [element for element in clean_formula if '+' in element]
    chemical_groups = [element for element in clean_formula if any([moiety in element for moiety in chemical_groups])]
    remainder = sorted([element for element in clean_formula if not any([element in adducts, element in customs, element in chemical_groups])])
    clean_formula = customs + remainder + chemical_groups + adducts
    return '_'.join(clean_formula)


def main():
    # Main function
    parser = argparse.ArgumentParser(description='Compare the output of mass_pattern_finder with a reference')
    parser.add_argument('-input', help='Path to the input file.', type=str)
    parser.add_argument('-reference', help='Path to reference file.', type=str)
    parser.add_argument('-time_tolerance', help='Tolerance in minutes for peaks to be comparable. Default = 0.05.', type=float, default=0.05)
    parser.add_argument('-mass_tolerance', help='Mass tolerance for masses to be equal. Default = 5 ppm.', default = 5e-6, type=float)
    parser.add_argument('-plot_time_range',help='Time range for plotting', default=0-15, type=str)
    parser.add_argument('-plot_mass_range',help='Mass range to use for plotting', default='0-1000', type=str)
    parser.add_argument('-min_intensity', help='Minimum intensity for unique plots to be reported. Default=1e4', type=float, default=1e4)
    parser.add_argument('-overwrite', help='If overwrite is True, the data saved from previous runs will be overwritten.', action='store_true')
    parser.add_argument('-output_dir', help='Optional other output directory.', type=str, default=None)
    args = parser.parse_args()

    if os.path.isdir(args.input):
        input_files = [os.path.join(args.input, file) for file in os.listdir(args.input) if '_analyzed.txt' in file.lower()]
    else:
        input_files = [os.path.join(os.path.dirname(args.input), file) for file in os.listdir(os.path.dirname(args.input)) if all(['_analyzed.txt' in file.lower(), os.path.basename(args.input) in file])]
    if os.path.isdir(args.reference):
        reference_files = [os.path.join(args.reference, file) for file in os.listdir(args.reference) if '_analyzed.txt' in file.lower()]
    else:
        reference_files = [os.path.join(os.path.dirname(args.reference), file) for file in os.listdir(os.path.dirname(args.reference)) if all(['_analyzed.txt' in file.lower(), os.path.basename(args.reference) in file])]

    sys.stdout.write(f"Obtaining peaks from analyzed files.\n{''.join(['=' for i in range(20)])}\n")
    input_peaks = parse_files(input_files, args, sample=True)
    reference_peaks = parse_files(reference_files, args, sample=False)

    if len(input_files) > 1:
        input_peaks = filter_peaks_for_identical_masses(input_peaks, args)
    if len(reference_files)>1:
        reference_peaks = filter_peaks_for_identical_masses(reference_peaks, args)

    sys.stdout.write(f"Comparing peaks.\n{''.join(['=' for i in range(20)])}\n")
    unique_peaks, difference_peaks = compare_peaks(input_peaks, reference_peaks, args, 'compare')

    if args.output_dir is not None and not os.path.isdir(args.output_dir):
        sys.stdout.write(f"ERROR: Could not find output directory {args.output_dir}. Saving in {os.path.dirname(args.input)} without overwriting.\n")
        args.output_dir = os.path.dirname(args.input)
        args.overwrite = False
    if args.output_dir is None:
        args.output_dir = os.path.dirname(args.input)

    output_file = f"{args.output_dir}/{os.path.splitext(os.path.basename(args.input))[0]}_ref_{os.path.splitext(os.path.basename(args.reference))[0]}.txt"
    with open(output_file, 'w') as file:
        file.write(f"Input file: {args.input}\nReference file: {args.reference}\n\nFound {len(unique_peaks)} unique peaks.\n\n")
        file.write(f"Time tolerance: {args.time_tolerance}\nMass tolerance: {args.mass_tolerance*1e6} ppm\n\n")
        file.write(''.join(['=' for _ in range(20)])+ '\n')
        for rt in unique_peaks:
            if unique_peaks[rt] == []:
                continue
            file.write(f"Unique peak at {rt}:\n")
            for mass in unique_peaks[rt]:
                file.write(f"\tMass:{mass[0]}\tIntensity: {mass[1]}\n\t\t")
                file.write('\n\t\t'.join([formula for formula in mass[2]])+'\n')
    sys.stdout.write(f"Saved unique peaks to {output_file}.\n")

    time_range = [float(time) for time in args.plot_time_range.split('-')]
    mass_range = [float(mass) for mass in args.plot_mass_range.split('-')]
    plot_filename = f"{args.output_dir}/{os.path.splitext(os.path.basename(args.input))[0]}_ref_{os.path.splitext(os.path.basename(args.reference))[0]}.png"
    plot_results_in_2D(unique_peaks, plot_filename, time_range, mass_range, args)

if __name__=='__main__':
    main()
