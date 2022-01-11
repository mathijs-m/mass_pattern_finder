# -*- coding: utf-8 -*-
"""
Script to compare output of the mass_pattern_finder with a reference file
@author: mmabesoone
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

def plot_results_in_2D(unique_peaks, plot_filename, time_range, mass_range, overwrite):
    # Plot the analyzed spectra in a single graph
    plot_list = list()
    for rt in unique_peaks:
        spectrum = unique_peaks[rt]
        for peak in spectrum[1:]:
            plot_list.append([peak[0], rt, log10(peak[1])])
    plot_list = pd.DataFrame(plot_list, columns=['mass', 'time', 'intensity'])
    plot_list = plot_list.sort_values(by='intensity', ascending='True', ignore_index=True)

    # Define the plot
    zrange = [power for power in range(floor(min(plot_list['intensity'])), ceil(max(plot_list['intensity'])))]
    plt.figure(figsize=(12, 12))
    plt.title('Unique masses\nInput: ' + '\nRef: '.join(os.path.basename(plot_filename).split('_ref_')), fontsize=16)
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
    sys.stdout.write(str(overwrite))
    if not overwrite and any([os.path.isfile(plot_filename+'.svg'), os.path.isfile(plot_filename+'.png')]):
        suffix = 1
        while any([os.path.isfile(plot_filename+'_'+str(suffix)+'.svg'), os.path.isfile(plot_filename+'_'+str(suffix)+'.png')]):
            suffix += 1
        plot_filename = f"{plot_filename}_{suffix}"
    plt.savefig(plot_filename+'.svg', transparent=True, dpi=300)
    plt.savefig(plot_filename+'.png', transparent=True, dpi=300)
    sys.stdout.write(f"Saved plot:\n\t{plot_filename}.svg\n\t{plot_filename}.svg\n")



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

    args = parser.parse_args()
    args.overwrite = bool(args.overwrite)
    input_peaks = obtain_peaks_from_file(args.input)
    reference_peaks = obtain_peaks_from_file(args.reference)
    
    unique_peaks = defaultdict(list)
    difference_peaks = defaultdict(list)

    num_reference_peaks = len(reference_peaks)
    elution_reference_peaks = [elution for elution in reference_peaks]
    ref_iterator = 0

    for elution_time in input_peaks:
        # Check if the peaks at all elution times in the sample match to a reference
        while elution_reference_peaks[ref_iterator]-elution_time < -1*args.time_tolerance:
            ref_iterator += 1
        if elution_reference_peaks[ref_iterator]-elution_time > args.time_tolerance:
            unique_peaks[elution_time] = input_peaks[elution_time]
            continue
        elif abs(elution_reference_peaks[ref_iterator]-elution_time) < args.time_tolerance:
           input_masses = input_peaks[elution_time]
           reference_masses = [mass for mass in reference_peaks[elution_reference_peaks[ref_iterator]]]
           for input_mass in input_masses:
               for reference_mass in reference_masses:
                   if (input_mass[0]-reference_mass[0]) < args.mass_tolerance*input_mass[0]:
                       difference_peaks[elution_time].append([input_mass[0], input_mass[1]-reference_mass[1], input_mass[2]])
               if not any([abs(input_mass[0]-reference_mass[0])<args.mass_tolerance*input_mass[0] for reference_mass in reference_masses]) and input_mass[1] > args.min_intensity:
                   unique_peaks[elution_time].append(input_mass)

    output_file = f"{os.path.splitext(args.input)[0]}_ref_{os.path.splitext(os.path.basename(args.reference))[0]}.txt"
    with open(output_file, 'w') as file:
        file.write(f"Input file: {args.input}\nReference file: {args.reference}\n\nFound {len(unique_peaks)} unique peaks.\n\n")
        for rt in unique_peaks:
            file.write(f"Unique peak at {rt}:\n")
            for mass in unique_peaks[rt]:
                file.write(f"\tMass:{mass[0]}\tIntensity: {mass[1]}\n\t\t")
                file.write('\n\t\t'.join([formula for formula in mass[2]])+'\n')
    sys.stdout.write(f"Saved unique peaks to {output_file}.\n")
    
    time_range = [float(time) for time in args.plot_time_range.split('-')]
    mass_range = [float(mass) for mass in args.plot_mass_range.split('-')]
    plot_filename = f"{os.path.splitext(args.input)[0]}_ref_{os.path.splitext(os.path.basename(args.reference))[0]}.png"
    plot_results_in_2D(unique_peaks, plot_filename, time_range, mass_range, args.overwrite)        
        
if __name__=='__main__':
    main()
