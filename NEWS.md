# 0.4.5

## Bug fixes

- Fixed both core and busco trees not being including in `plot_tree`

# 0.4.4

## Bug fixes

- Fixed invalid `ggtree` dependency error

# 0.4.3

## Bug fixes

- Fixed NA values for ANI in `estimated_ani_match_table`

# 0.4.2

## Improvements

- Replaced tree plotting functionality with `heattree`

# 0.4.1

## New features

- Added `print_outdir_schema()` function for displaying output directory schema information
- Added `verbose` and `print_all` options to `print_outdir_schema()`
- Added `must_exist` parameter to `find_ps_data()` function

## Improvements

- Replaced the `subset` option with `max_size` and `prefer_contextual` parameters for better control
- Improved error messages in `find_ps_paths()` function
- Enhanced output metadata handling with name and description fallbacks

