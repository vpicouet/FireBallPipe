#! /usr/bin/env python2
from __future__ import print_function
from __future__ import division


import os, sys
import numpy as np
import glob
from astropy.io import fits


def ChangeNameGuiderImages(path):
    print(path)
    files = glob.glob(path + "*.fits")
    print(files)
    print("\nReading fits files...")
    for i, file in enumerate(files):
        with fits.open(file) as f:
            name = f.filename()
            rot = "%+04d" % (f[0].header["ROTENC"])
            lin = ""  #'%04d'%(100*f[0].header['LINAENC'])
            date = f[0].header["DATE"].replace(":", "-")
            os.rename(
                name, f.filename()[:-5] + "_pa" + rot + "" + lin + "_" + date + ".fits"
            )
    return 7971390


# os.rename('stack8103901.fits', 'stack8103901_test.fits')


def revert_names(path):
    print(path)
    files = glob.glob(path + "*.fits")
    print(files)
    print("\nReading fits files...")
    for i, file in enumerate(files):
        with fits.open(file) as f:
            name = f.filename()
            if len(os.path.basename(name)) == 17:
                # print ('Names have already a nominal shape')
                raise Exception("Names have already a nominal shape")

            newname = (
                os.path.dirname(name)
                + "/"
                + os.path.basename(name)[:12]
                + os.path.basename(name)[-5:]
            )
            print(name, "\n", newname)  # + name[:12] + name[-5:])
            os.rename(name, newname)
    return


if __name__ == "__main__":
    print(
        """
        *************************************
        *     Program change name guider images    *
        *************************************
        """
    )

    import argparse

    parser = argparse.ArgumentParser(prog="ChangeNames")
    parser.add_argument(
        "-p",
        "--path",
        default=os.getcwd(),  # type='string',
        help="""The path of the path to change the guider images""",
    )
    parser.add_argument(
        "-r",
        "--reverse",
        action="store_true",  # type='string',
        help="""put the old names""",
    )

    args = parser.parse_args()
    print(args.reverse)

    if args.reverse:
        revert_names(args.path)
    else:
        ChangeNameGuiderImages(args.path)


#    print(dedent("""
#        *************************************
#        *     Program process_catalog.py    *
#        *************************************
#        """))
#
#    density_parser = density.parser()
#
#    density_parser.set_defaults(NDnetField='dtfe')
#
#    parser = argparse.ArgumentParser(prog="process_catalog", parents=[density_parser],)
#
#    parser.add_argument('-s', '--nsigma', dest='nsig', default=3.0, type=float,
#                        help='persistense nsigma threshold, (default: %(default)s)')
#
#    parser.add_argument('-c', '--cut', dest='cutp', type=float,
#                        help='persistense threshold, (nsigma threshold used by default)')
#
#    parser.add_argument('-k', '--smooth_skl', dest='nskel_smooth', default=1, type=int,
#                        help='smooth output skeleton and walls (default: %(default)s)')
#
#    parser.add_argument('-o', '--out', dest='outnamebase',
#                        help='output name prefix (eventually including the path)')
#
#    parser.add_argument('--skl', dest='sklfilename',
#            help='''use an already computed skeleton stored in the file sklfilename.
#                    The skeleton should have been computed using this command
#                    & the name is the one of the NDskl ascii file
#                    (without the extension .a.NDskl)''')
#
#    parser.add_argument('--training', action='store_true',
#                        help='compute the skeleton only and plot the persistence diagrams')
#
#    parser.add_argument('-f', '--full', action='store_true',
#                            help='same than --fof --walls --patches')
#
#    parser.add_argument('--fof',  action='store_true',
#                        help='''compute groups arround maxima using a lower limit in density.
#                                Do not compute them if not present (default).''')
#
#    parser.add_argument('--fof_threshold', default=0., type=float,
#            help= '''set the lower limit in density for computing fofs. It is expressed as a
#            a ratio between the density of the highest saddle connected to the maxima and
#            the density of the maxima:threshold = (1-ratio) * saddle + ratio * peak
#            i.e. 0. (default) means density of the saddle, 1. means density of the maxima''')
#
#    parser.add_argument('--walls', action='store_true',
#                        help='''compute distance of galaxies to the nearest wall.
#                                Do not compute them if not present (default).''')
#
#    parser.add_argument('--patches', action='store_true',
#                        help='''compute voids and peaks patches.
#             (determine for each galaxy the index of the min and max whose patch it belongs to)
#             Do not compute them if not present (default).''')
#
#    args = parser.parse_args()
#
#
#    ## build default output name
#    catfilenamebase, ext = path.splitext(args.catfilename)
#    if not args.outnamebase:
#        args.outnamebase = catfilenamebase
#
#    if not args.sklfilename:
#        if args.weight: weight_str = '_weight'
#        else: weight_str = ''
#
#        if args.NDnetField == 'dtfe': pstfix = '_dtfe' + weight_str
#        elif args.NDnetField == 'log_dtfe': pstfix = '_log_dtfe' + weight_str
#        elif args.NDnetField == 'haloes':
#            pstfix = '_haloes' + weight_str
#            pstfix += '_SZ{0:.3g}'.format(args.haloes_size)
#        elif args.NDnetField == 'log_haloes':
#            pstfix = '_log_haloes' + weight_str
#            pstfix += '_SZ{0:.3g}'.format(args.haloes_size)
#        elif args.NDnetField == 'gkernel':
#            pstfix = '_gkernel' + weight_str
#            pstfix += '_SZ{0:.3g}'.format(args.gkernel_sigma)
#        elif args.NDnetField == 'log_gkernel':
#            pstfix = '_log_gkernel' + weight_str
#            pstfix += '_SZ{0:.3g}'.format(args.gkernel_sigma)
#        elif args.NDnetField == 'sphere':
#            pstfix = '_sphere' + weight_str
#            pstfix += '_SZ{0:.3g}'.format(args.sphere_radius)
#        elif args.NDnetField == 'log_sphere':
#            pstfix = '_log_sphere' + weight_str
#            pstfix += '_SZ{0:.3g}'.format(args.sphere_radius)
#        elif args.NDnetField == 'neighbors':
#            pstfix = '_neighbors' + weight_str
#            pstfix += '_NB{0:d}'.format(args.neighbors_nb)
#        elif args.NDnetField == 'log_neighbors':
#            pstfix = '_log_neighbors' + weight_str
#            pstfix += '_NB{0:d}'.format(args.neighbors_nb)
#
#        if args.cutp:
#            pstfix += '_SD{0:d}_C{1:.3g}_K{2:d}'.format(args.ndens_smooth, args.cutp, args.nskel_smooth)
#        else:
#            pstfix += '_SD{0:d}_S{1:g}_K{2:d}'.format(args.ndens_smooth, args.nsig, args.nskel_smooth)
#    else:
#        pstfix = '_skel_' + path.basename(args.sklfilename)
#    args.outnamebase += pstfix
#
#    # handle full optional argument
#    if args.full:
#        args.fof = True
#        args.walls = True
#        args.patches = True
#    del args.full
#
#    print("Called with following argument & options:\n", args, "\n")
#
#    # run it
#    cat, delaunay_cat = process_catalog(**vars(args))
#
#    ## write final vtk & ascii catalog (delaunay grid + added data)
#    ############################################################
#    delaunay_cat.write_vtu(args.outnamebase + ".vtu")
#
#    cat.write(args.outnamebase + ext, overwrite=True)  # write with everything
#
