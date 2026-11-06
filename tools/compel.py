#!/usr/bin/python

"""
Compel 0.01 2026-11-05
  Fit the minimum-volume ellipsoid enclosing well-measured reflections in an MTZ file.

  Reflections used for the fit: valid FP and SIGFP with FP/SIGFP > --min-snr (default 2.0),
  expanded by the space-group operators (including Friedel's law). Reports resolution limits
  along a*, b*, c* and the fraction of reflections inside the ellipsoid that have valid
  FP and SIGFP. With --trim, writes a new MTZ file without the reflections that lie outside
  the ellipsoid and have no valid FP/SIGFP.

  Written by Robbie Joosten with help of AI
  E-mail:    r.joosten@nki.nl

  If you publish results (directly or indirectly) obtained by using
  Compel: Please, refer to (one of) these references:
  - Robbie P. Joosten, Jean Salzemann, Vincent Bloch, Heinz Stockinger,
    Ann-Charlott Berglund, Christophe Blanchet, Erik Bongcam-Rudloff,
    Christophe Combet, Ana L. Da Costa, Gilbert Deleage, Matteo
    Diarena, Roberto Fabbretti, Geraldine Fettahi, Volker Flegel,
    Andreas Gisel, Vinod Kasam, Timo Kervinen, Eija Korpelainen, Kimmo
    Mattila, Marco Pagni, Matthieu Reichstadt, Vincent Breton, Ian J.
    Tickle, Gert Vriend: "PDB_REDO: automated re-refinement of X-ray
    structure models in the PDB" J. Appl. Cryst., 42, p. 376-384 (2009)
"""

import argparse

import gemmi
import numpy as np
from scipy.spatial import ConvexHull

__version__ = '0.01'


def fit_ellipsoid(points, eps=1e-4, max_iter=100000):
    """Return the matrix A of the minimum-volume ellipsoid {h : h^T A h <= 1} enclosing points.

    Uses Khachiyan's algorithm, which iteratively reweights the points until the
    ellipsoid is optimal to within the tolerance eps. Only the convex hull vertices
    matter for the enclosing ellipsoid, so the hull is computed first to speed things up.
    """
    hull = ConvexHull(points)
    P = points[hull.vertices]
    n, dim = P.shape
    # Homogeneous coordinates [h, 1] used by the algorithm
    Q = np.hstack([P, np.ones((n, 1))]).T
    # Weights on each hull vertex, starting uniform
    u = np.ones(n) / n
    for _ in range(max_iter):
        X = (Q * u) @ Q.T
        # Leverage of each point; the optimum has max leverage equal to dim + 1
        leverage = np.einsum('ij,ji->i', Q.T @ np.linalg.inv(X), Q)
        j = np.argmax(leverage)
        if leverage[j] <= (1 + eps) * (dim + 1):
            break
        # Move weight towards the point that is furthest outside the current ellipsoid
        step = (leverage[j] - dim - 1) / ((dim + 1) * (leverage[j] - 1))
        u = (1 - step) * u
        u[j] += step
    # Centre and shape matrix from the final weights; the centre is ~0 by Friedel symmetry
    centre = P.T @ u
    A = np.linalg.inv(P.T @ np.diag(u) @ P - np.outer(centre, centre)) / dim
    return A


def main():
    parser = argparse.ArgumentParser(description='Fit an ellipsoid enclosing well-measured reflections.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('mtz', help='input MTZ file')
    parser.add_argument('--trim', metavar='OUTFILE',
                        help='write a new MTZ file without reflections outside the ellipsoid that lack valid FP and SIGFP')
    parser.add_argument('--min-snr', type=float, default=2.0,
                        help='minimum FP/SIGFP for a reflection to be used in the fit (default 2.0)')
    args = parser.parse_args()

    # Read the MTZ file and pull out the Miller indices and the FP/SIGFP columns
    mtz = gemmi.read_mtz_file(args.mtz)
    data = np.array(mtz, copy=True)
    idx = {label: i for i, label in enumerate(mtz.column_labels())}
    hkl = data[:, [idx['H'], idx['K'], idx['L']]].astype(float)
    fp = data[:, idx['FP']]
    sig = data[:, idx['SIGFP']]

    # A reflection is valid if both FP and SIGFP are present; strong if also above the S/N threshold
    valid = ~np.isnan(fp) & ~np.isnan(sig)
    with np.errstate(invalid='ignore', divide='ignore'):
        strong = valid & (fp / sig > args.min_snr)

    # Symmetry expansion: apply every space-group rotation, then each combined with Friedel's law (-I),
    # so the fitted set of reflections is centrosymmetric and closed under the point group
    rotations = [np.array(op.rot, dtype=float) / gemmi.Op.DEN for op in mtz.spacegroup.operations()]
    rotations += [-R for R in rotations]
    fit_hkl = hkl[strong]
    expanded = np.vstack([fit_hkl @ R for R in rotations])
    # Remove duplicates created by symmetry-equivalent reflections
    expanded = np.unique(np.round(expanded).astype(int), axis=0).astype(float)

    A = fit_ellipsoid(expanded)
    # Rescale so the largest value of h^T A h over the fitted reflections is exactly 1,
    # guaranteeing that every fitted reflection lies inside the ellipsoid
    A = A / np.einsum('ij,jk,ik->i', expanded, A, expanded).max()

    # Report the resolution limit along each reciprocal axis: the ellipsoid meets the h-axis
    # at h = 1/sqrt(A_hh), and d = 1/(h * |a*|)
    rec = mtz.cell.reciprocal()
    print(f'File: {args.mtz}')
    print(f'Space group: {mtz.spacegroup.hm}')
    print(f'Reflections used for fit (FP/SIGFP > {args.min_snr}): {int(strong.sum())}')
    print(f'Symmetry-expanded reflections: {len(expanded)}')
    print()
    print('Resolution along reciprocal axes:')
    for i, (name, length) in enumerate(zip(['a*', 'b*', 'c*'], [rec.a, rec.b, rec.c])):
        h_max = 1 / np.sqrt(A[i, i])
        d_min = 1 / (h_max * length)
        print(f'  {name}: index limit {h_max:.2f}, d_min {d_min:.3f} A')

    # Test every reflection in the file against the ellipsoid (inside means h^T A h <= 1)
    all_hkl = hkl
    inside = np.einsum('ij,jk,ik->i', all_hkl, A, all_hkl) <= 1
    n_inside = int(inside.sum())
    n_valid = int((inside & valid).sum())
    print()
    print(f'Reflections inside ellipsoid: {n_inside}')
    print(f'With valid FP and SIGFP: {n_valid}')
    print(f'Fraction with valid FP and SIGFP: {100 * n_valid / n_inside:.1f}%')

    if args.trim:
        # Keep reflections inside the ellipsoid, and any reflection that has valid values
        keep = inside | valid
        mtz.set_data(data[keep])
        mtz.write_to_file(args.trim)
        print()
        print(f'Removed {int((~keep).sum())} reflections outside the ellipsoid without valid FP and SIGFP')
        print(f'Wrote {int(keep.sum())} reflections to {args.trim}')


if __name__ == '__main__':
    main()
