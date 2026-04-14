# Peak Finder

Desktop tooling for inspecting mass spectrometry peaks, comparing manual and automated peak bounds, and exporting review results to CSV.

The main application is `peak_inspector.py`, a Tkinter GUI that loads a two-column spectrum file (`m/z`, `intensity`), lets you step through target peaks, and compares several automatic boundary-selection methods against manual picks.

## Project layout

```
peak_finder/
├── peak_inspector.py          # main interactive desktop application
├── auto_selection.py          # core automatic peak-selection algorithms
├── data/
│   ├── reference/
│   │   ├── isotopes.json      # parsed isotope database (generated from NIST source)
│   │   └── isotopic_nist_info.txt  # raw NIST isotope source data
│   └── samples/
│       ├── neodynium.csv      # default example spectrum (loaded on startup)
│       └── excess noiseavg scan.csv  # additional sample spectrum
├── examples/
│   ├── test.csv               # example exported summary output
│   └── test_methods.csv       # example exported per-method output
├── tools/
│   ├── txt_to_json.py         # converts isotopic_nist_info.txt → isotopes.json
│   └── peak_finder.py         # standalone derivative plotting script
└── logs/                      # runtime logs, one file per session (git-ignored)
```

## What it does

- Loads a spectrum from CSV and analyzes one or more target `m/z` values.
- Applies optional Savitzky-Golay smoothing before automated peak detection.
- Computes automatic bounds using derivative, threshold, and minima-based methods.
- Supports manual left / peak / right assignment for review and comparison.
- Exports a summary CSV plus a per-method CSV for downstream analysis.
- Matches reviewed peaks against isotope combinations from `data/reference/isotopes.json`.
- Writes runtime logs to the `logs/` directory.

## Requirements

Python 3.10+ is a reasonable baseline.

Install the packages used by the project:

```powershell
pip install numpy pandas matplotlib scipy tksheet
```

Notes:

- `tkinter` is part of the Python standard library on most Windows Python installs.
- The GUI uses the Matplotlib `TkAgg` backend.

## Quick start

From the project folder:

```powershell
python peak_inspector.py
```

If you are using the local virtual environment already present in this workspace:

```powershell
.\venv\Scripts\Activate.ps1
python peak_inspector.py
```

On startup, the app will try to:

- load `data/samples/neodynium.csv`
- add the built-in default target list
- preload isotope data from `data/reference/isotopes.json`

## Input formats

### Spectrum CSV

The main spectrum file should contain two columns:

1. `m/z`
2. intensity

The loader accepts files with or without a header row and uses the first two columns.

### Target list

You can load a `.txt` or `.csv` target list from the File menu.

- One target `m/z` per line is expected.
- Comment lines beginning with `#` are ignored.
- For CSV-style lines, the first comma-separated value is used.

Example:

```text
# neodymium targets
157.901
158.903
159.904
```

## GUI workflow

Typical review flow:

1. Launch `peak_inspector.py`.
2. Load a spectrum CSV if you do not want the default example file.
3. Load targets from file or use the built-in defaults.
4. Navigate through targets and inspect the automatic peak/bound selections.
5. Assign manual left, peak, and right positions when needed.
6. Export the results from the File menu.

The application can compare several automatic methods:

- derivative bounds
- threshold bounds
- minima bounds
- combined bounds

It also exposes analysis preferences for smoothing, window sizes, and active method selection through the Analysis menu.

Spreadsheet-style table selection is available throughout the main grid views.

- Click cells, row headers, or column headers to select like a spreadsheet.
- Copy selections from the table views directly as TSV with `Ctrl+C`.
- Paste TSV into the Targets table to add masses from the clipboard.
- Paste isotope labels into the selected-isotopes table in the isotope selector to add matching labels from the loaded isotope database.

## Keyboard shortcuts

Inside the main window:

- Right arrow: next target
- Left arrow: previous target
- `Ctrl+C`: copy the current table selection as TSV
- `Ctrl+V`: paste TSV into the active supported table view
- `l`: assign manual left bound
- `p`: assign manual peak center
- `r`: assign manual right bound
- `c`: clear manual points for the current target

## Exported results

When you export results, the app writes two CSV files:

- the summary file you choose, for example `results.csv`
- a second file beside it named `results_methods.csv`

The exported data includes, among other fields:

- target `m/z`
- automatic and manual peak positions
- automatic and manual bounds
- bound widths
- integrated intensity within each bound set
- percent error versus manual selection
- best isotope label match, with error in `mmu` and `ppm`

See `examples/test.csv` and `examples/test_methods.csv` for the expected output format.

## Isotope database

The GUI uses `data/reference/isotopes.json` to support isotope labeling and combination matching.

To regenerate that JSON from the bundled NIST-style text file:

```powershell
python tools/txt_to_json.py
```

This reads `data/reference/isotopic_nist_info.txt` and writes `data/reference/isotopes.json`.

## Logs

Application logs are written to `logs/` with timestamped filenames such as:

```text
logs/peak_inspector_YYYYMMDD_HHMMSS.log
```

These logs are useful if the GUI hits an unexpected exception. The `logs/` directory is git-ignored.

## Tools

`tools/peak_finder.py` is a small standalone plotting script for exploring a narrow region of the neodymium spectrum and visualizing the raw signal plus first and second derivatives. It is separate from the main GUI workflow.

`tools/txt_to_json.py` converts the NIST isotope text file into the JSON database used by the GUI. Run it from the project root whenever `data/reference/isotopic_nist_info.txt` is updated.

## Known assumptions

- The spectrum loader expects numeric data in the first two columns.
- Automatic analysis works on a local window around each target `m/z`.
- If `data/reference/isotopes.json` is missing or invalid, isotope matching is simply unavailable rather than fatal.
- Run `peak_inspector.py` from the project root so relative data paths resolve correctly.
