# msPepSeq

## Program for Peptide Sequencing based on b- and y-ion Fragmentation

This program performs peptide sequencing from mass spectrometry data using b- and y-ion fragmentation. It assigns possible peptide sequences to observed m/z peaks and reports candidate sequences that best match the input spectrum.

## Usage

```bash
python msPepSeq.py <path-to-ms-spectra> [--output_to_file] [--output_assignments]
```

- `<path-to-ms-spectra>`: Path to a text file containing one m/z value per line (e.g. `testData/b_y_spectrum.txt`).
- `--output_to_file`: (optional) Outputs the candidate sequences to a file in the `out/` directory instead of the console.
- `--output_assignments`: (optional) Saves all calculated assignments per peak for b- and y-ions to separate files in the `out/` directory.

Example:
```bash
python msPepSeq.py testData/b_y_spectrum.txt --output_to_file --output_assignments
```

## Notes

- The program is based exclusively on b- and y-ion fragmentation (no other fragment types are considered).
- Monoisotopic masses for amino acids are taken from my systemsbiology lecture.
- The mass tolerance is set to 0.055 Da.
- Output is either printed to the console or written to the `out/` directory (created automatically if needed).
- Candidates are ranked by similarity to the measured spectrum (based on Euclidean distance).

---
See the source code for further details on the algorithms and implementation.


