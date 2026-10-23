#!/usr/bin/env python3
"""
Convert an EAM vertical-coordinate netCDF file (containing hyai, hybi, hyam,
hybm) into the pair of standalone HOMME ascii vertical coordinate files
(e.g. acme-72i.ascii / acme-72m.ascii).

Usage:
    ./nc2ascii_vcoord.py input.nc output_prefix

Produces:
    output_prefix + "i.ascii"   (interface coefficients: hyai, hybi)
    output_prefix + "m.ascii"   (midpoint coefficients:  hyam, hybm)
"""
import sys
import argparse
import netCDF4 as nc


def write_section(f, values, label):
    n = len(values)
    f.write("{:>7d}   ! {}\n".format(n, label))
    for v in values:
        f.write("{:.15g}\n".format(float(v)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input EAM vertical coordinate netCDF file")
    parser.add_argument("outprefix", help="output file prefix (e.g. acme-72)")
    args = parser.parse_args()

    ds = nc.Dataset(args.infile)
    try:
        hyai = ds.variables["hyai"][:]
        hybi = ds.variables["hybi"][:]
        hyam = ds.variables["hyam"][:]
        hybm = ds.variables["hybm"][:]
    except KeyError as e:
        sys.exit(f"Error: missing expected variable {e} in {args.infile}")

    nlevi = len(hyai)
    nlevm = len(hyam)
    if len(hybi) != nlevi or len(hybm) != nlevm or nlevi != nlevm + 1:
        sys.exit(
            f"Error: inconsistent dimensions "
            f"(hyai={len(hyai)}, hybi={len(hybi)}, hyam={len(hyam)}, hybm={len(hybm)})"
        )

    iname = args.outprefix + "i.ascii"
    mname = args.outprefix + "m.ascii"

    with open(iname, "w") as f:
        write_section(f, hyai, "hyai")
        write_section(f, hybi, "hybi")

    with open(mname, "w") as f:
        write_section(f, hyam, "hyam")
        write_section(f, hybm, "hybm")

    print(f"Wrote {iname} ({nlevi} interfaces)")
    print(f"Wrote {mname} ({nlevm} levels)")


if __name__ == "__main__":
    main()
