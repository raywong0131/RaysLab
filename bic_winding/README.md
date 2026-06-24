# BIC winding number calculator

This small CLI computes the winding number of a BIC polarization singularity
from tabular `kx, ky, cx, cy` data.

Default input convention:

- column B: `kx`
- column C: `ky`
- column D: `cx`
- column E: `cy`
- no header row

For grid data, `--path auto` extracts the outer rectangular boundary and uses
that as the closed loop. For data already sampled along a loop, use
`--path input`. For unordered loop samples, use `--path angle`.

## Run on the provided example

```powershell
& "C:\Users\wangs\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" `
  "D:\CodeXprojects\bic_winding\bic_winding.py" `
  "D:\Claude Code\PolarSingularity\K(-0.1-0.1)polarization for b245eta0.96sca1.156.xlsx"
```

Expected result for the example workbook with the default `auto` settings:

```text
winding number   : -1
nearest integer  : -1
```

## Field modes

- `auto`: use `complex-axis` when `cx, cy` have significant imaginary parts,
  otherwise use `real`.
- `real`: compute `arg(real(cx) + i real(cy))`.
- `imag`: compute `arg(imag(cx) + i imag(cy))`.
- `complex-axis`: treat `cx, cy` as a complex Jones vector, remove the local
  common phase, and compute the winding of the recovered real polarization
  axis.
- `stokes-axis`: compute the polarization-axis angle from Stokes-like
  quantities. This is useful for checking complex Jones-vector data.

The sign of the winding depends on the loop direction. Add `--reverse` if your
sampling convention is clockwise and you want the opposite orientation.

`--plot output.png` can create a diagnostic plot when `matplotlib` is available.

Use `--scale-x`, `--scale-y`, and `--half-width` to manually rescale the vector
components and choose a square loop before computing the winding number, for
example `--scale-y 10 --half-width 0.02`.

The same values can be edited directly near the top of `bic_winding.py`:

```python
# half_width = 0.02
# scale_x = 1.0
# scale_y = 10.0
```

If a line stays commented out, the default is used. The default scales are
`scale_x = 1` and `scale_y = 1`; the default `half_width` behavior is automatic.

For rectangular grid data, the program also performs a shrink scan: it starts
from the largest square loop around the selected center, reduces `half_width`
by the available grid spacing, recomputes the winding number on each loop, and
reports the first `half_width` where the winding changes.
