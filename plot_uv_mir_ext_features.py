import argparse

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec

from astropy.table import Table
import astropy.units as u

from dust_extinction.averages import GCC09_MWAvg
from measure_extinction.extdata import ExtData

if __name__ == '__main__':

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    # setup plot

    fontsize = 18

    font = {'size': fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('ytick.minor', width=2)

    fig, ax = plt.subplots(figsize=(18, 7))
    gs = gridspec.GridSpec(1, 11)
    ax = []
    ax.append(plt.subplot(gs[0, 0:4]))
    ax.append(plt.subplot(gs[0, 4:7]))
    ax.append(plt.subplot(gs[0, 7:11]))

    # plot the MW average
    ext_model = GCC09_MWAvg()
    ax[0].plot(1.0/ext_model.obsdata_x_fuse, ext_model.obsdata_axav_fuse, 'b-')
    ax[0].plot(1.0/ext_model.obsdata_x_iue, ext_model.obsdata_axav_iue, 'b-')
    ax[0].text(0.22, 3.25, r'$2175$ $\AA$', horizontalalignment='center')
    ax[0].text(0.12, 4.5, r'FUV rise', horizontalalignment='center',
               rotation=-55.)
    wave_bands = [0.11, 0.15, 0.19, 0.2175, 0.25]
    frac_width = 0.05
    for cwave in wave_bands:
        bx = [cwave-frac_width*cwave, cwave-frac_width*cwave,
              cwave+frac_width*cwave, cwave+frac_width*cwave]
        by = [0.0, 1.0, 1.0, 0.0]
        ax[0].plot(bx, by, 'k-')
    ax[0].text(0.175, 1.4, 'desired spectral range',
               horizontalalignment='center', fontsize=fontsize*0.75)
    ax[0].errorbar([0.2], [1.3], xerr=[0.1], fmt='o',
                   markersize=0.0001, capsize=10,
                   elinewidth=2)
    ax[0].text(0.175, 0.5, 'desired imaging bands',
               horizontalalignment='center', fontsize=fontsize*0.75,
               bbox=dict(facecolor='white'))

    ax[0].set_yscale('linear')
    ax[0].set_xscale('log')
    ax[0].set_ylim(0.0, 6.5)
    ax[0].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[0].set_ylabel(r'$A(\lambda)/A(V)$')
    ax[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax[0].get_xaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax[0].set_xticks([0.1, 0.2, 0.3])
    ax[0].set_xticks([0.15], minor=True)

    # read in and plot the 3.4 um feature
    a = Table.read('vicyg12_34um_fromplot.txt', format='ascii.basic')
    a['wave'].unit = '1/cm'
    freq = a['wave'].quantity
    waves = freq.to(u.micron, equivalencies=u.spectral())
    tau = a['tau']/10.2 + 0.075
    ax[1].plot(waves, tau, 'r-')

    ax[1].text(3.4, 0.097, r'3.4 $\mu$m', horizontalalignment='center')
    wave_bands = [3.2, 3.4, 3.6]
    frac_width = 0.02
    for cwave in wave_bands:
        bx = [cwave-frac_width*cwave, cwave-frac_width*cwave,
              cwave+frac_width*cwave, cwave+frac_width*cwave]
        by = [0.06, 0.066, 0.066, 0.06]
        ax[1].plot(bx, by, 'k-')
    ax[1].text(3.4, 0.0685, 'desired spectral range',
               horizontalalignment='center', fontsize=fontsize*0.75)
    ax[1].errorbar([3.4], [0.068], xerr=[0.3], fmt='o',
                   markersize=0.0001, capsize=10,
                   elinewidth=2)
    ax[1].text(3.4, 0.063, 'desired imaging bands',
               horizontalalignment='center', fontsize=fontsize*0.75,
               bbox=dict(facecolor='white'))

    ax[1].set_yscale('linear')
    ax[1].set_xscale('log')
    ax[1].set_xlim(3.0, 3.8)
    ax[1].set_ylim(0.06, 0.1)
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax[1].get_xaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax[1].set_xticks([3, 3.5])
    ax[1].set_xticks([3.25, 3.75], minor=True)

    # read in and plot the spitzer average extinction
    spitext = ExtData()
    spitext.read_ext_data('all_ext_trim_ave.fits')
    spitext.plot_ext(ax[2], color='r',
                     fontsize=fontsize,
                     legend_key='IRS',
                     legend_label='Average (this work)')
    wave_bands = [7.0, 10.0, 14.0]
    frac_width = 0.1
    for cwave in wave_bands:
        bx = [cwave-frac_width*cwave, cwave-frac_width*cwave,
              cwave+frac_width*cwave, cwave+frac_width*cwave]
        by = [-0.05, -0.025, -0.025, -0.05]
        ax[2].plot(bx, by, 'k-')
    ax[2].text(12.5, -0.018, 'desired spectral range',
               horizontalalignment='center', fontsize=fontsize*0.75)
    ax[2].errorbar([17.5], [-0.02], xerr=[12.5], fmt='o',
                   markersize=0.0001, capsize=10,
                   elinewidth=2)
    ax[2].text(9.75, -0.0375, 'desired imaging bands',
               horizontalalignment='center', fontsize=fontsize*0.75,
               bbox=dict(facecolor='white'))

    ax[2].text(10.0, 0.08, r'10 $\mu$m', horizontalalignment='center')
    ax[2].text(20.0, 0.045, r'20 $\mu$m', horizontalalignment='center')
    ax[2].set_yscale('linear')
    ax[2].set_xscale('log')
    ax[2].set_xlim(4.75, 32.0)
    ax[2].set_ylim(-0.05, 0.1)
    ax[2].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax[2].get_xaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax[2].set_xticks([10])
    ax[2].set_xticks([5, 6, 7, 8, 15, 20, 30], minor=True)

    plt.tight_layout(w_pad=0.5, h_pad=1.0)

    save_filename = 'uv_mir_ext_features'
    if args.png:
        fig.savefig(save_filename+'.png')
    elif args.pdf:
        fig.savefig(save_filename+'.pdf')
    else:
        plt.show()
