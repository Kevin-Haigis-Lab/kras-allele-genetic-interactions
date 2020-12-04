#!/bin/sh

for file in lib src tests munge 'paper/figures/make_figure_scripts'
do
	Rscript -e "styler::style_dir('$file')"
	# echo $file
done
