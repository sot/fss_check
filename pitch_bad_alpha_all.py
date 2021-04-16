import argparse
from pathlib import Path
from astropy.table import Table, vstack
from cxotime import CxoTime
import check_fss


def get_options():
    parser = argparse.ArgumentParser(description="Make long plot")
    parser.add_argument("--out",
                        default=".",
                        help="output directory, default is '.'")
    opt = parser.parse_args()
    return opt


def main():
    opt = get_options()
    outdir = opt.out
    if not Path(outdir).exists():
        Path(outdir).mkdir(parents=True)

    start = 2010
    stop = int(CxoTime().frac_year) + 1

    files = []
    for year in range(start, stop):
        fname = Path(outdir) / f'prim_dat_{year}.h5'
        if not Path(fname).exists() or year == stop - 1:
            dat = check_fss.get_fss_prim_data(start=f'{year}:001',
                                              stop=f'{year + 1}:001')
            Table(dat).write(fname, overwrite=True)
        files.append(fname)

    bigdat = vstack([Table.read(file) for file in files])
    check_fss.plot_pitches_any_kalman(bigdat,
                                      savefig=True,
                                      midrange=False,
                                      outdir=outdir)


if __name__ == '__main__':
    main()
