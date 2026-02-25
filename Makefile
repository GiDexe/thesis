.PHONY: all std no_std clean

# Default: All
all: output/preliminary_results.pdf

# --- NO_STD (N=32) ---
no_std: data/output_data/grid_search_ISO37001.RData \
        data/output_data/grid_search_ISO45001.RData

data/output_data/data_balanced_37001.rds \
data/output_data/data_balanced_45001.rds: cleaning_no_LOCF.r
	Rscript cleaning_no_LOCF.r

data/output_data/grid_search_ISO37001.RData \
data/output_data/grid_search_ISO45001.RData: grid_search_32.r \
	data/output_data/data_balanced_37001.rds \
	data/output_data/data_balanced_45001.rds
	Rscript grid_search_32.r

# --- STD (N=37) ---
std: data/output_data/grid_search_ISO37001_standardised.RData \
     data/output_data/grid_search_ISO45001_standardised.RData

data/output_data/grid_search_ISO37001_standardised.RData \
data/output_data/grid_search_ISO45001_standardised.RData: gird_search_standardised.r \
	data/output_data/data_balanced_37001.rds \
	data/output_data/data_balanced_45001.rds
	Rscript gird_search_standardised.r

# --- PRELIMINARY OUTPUT ---
output/preliminary_results.pdf: preliminary_results.qmd \
	data/output_data/grid_search_ISO37001.RData \
	data/output_data/grid_search_ISO45001.RData \
	data/output_data/grid_search_ISO37001_standardised.RData \
	data/output_data/grid_search_ISO45001_standardised.RData
	quarto render preliminary_results.qmd --output-dir output

# --- CLEAN ---
clean:
	rm -f output/tabs/*.tex output/tabs/*.pdf output/preliminary_results.pdf
	rm -f data/output_data/grid_search_*.RData