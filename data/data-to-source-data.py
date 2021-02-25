#!/bin/env python3

# A convoluted, yet safe way to create an archive of the data for submission.
# This script is not needed for the analysis and is not intended to be run by others.

import pathlib
from pathlib import Path
import shutil
from colorama import init, Fore, Back, Style
import pretty_errors
import getpass
import zipfile
import warnings

init(autoreset=True)

# If this is not being run by Joshua Cook, a User Warning is presented.
if getpass.getuser() != "jc604":
	warnings.warn("This script does not need to be run for the analysis.", UserWarning)

def info(s: str):
	print(Fore.BLUE + Style.DIM + s)
	return

def done():
	print(Fore.GREEN + Style.BRIGHT + "Done!")
	return

def error(s: str):
	print(Fore.RED + s)

DATA_DIR = Path("data")
TEMP_DIR = Path("/n/scratch3/users/j/jc604/Cook-et-al_2021_rawdata")

if not TEMP_DIR.exists():
	try:
		info("Creating temporary directory.")
		TEMP_DIR.mkdir()
	except FileExistsError as error:
		info("Temporary directory already exists.")
	except FileNotFoundError:
		info(f"Uh-oh. Not able to create the temporary directory as {TEMP_DIR}")
else:
	info("Temporary directory exists.")

directories_to_copy = [
	"annovar",
	"cancer-rates",
	"cbioportal",
	"depmap20Q1",
	"fisher-comutation",
	"kegg-pathways",
	"mutational-signatures",
	"tissue-gene-expression",
]


info("Copying data directories to a temporary destination.")
for d in directories_to_copy:
	src_dir = DATA_DIR / d
	dest_dir = TEMP_DIR / d
	if not src_dir.exists():
		info(f"... skipping '{d}' - does not exists.")
	elif dest_dir.exists():
		info(f"... skipping '{d}' - copy already exists.")
	else:
		try:
			info(f"... copying '{d}'")
			shutil.copytree(src_dir, dest_dir)
		except:
			error(f"Unable to copy directory '{d}'")
		

info("Zipping temporary directory.")
shutil.make_archive(base_name=DATA_DIR / TEMP_DIR.name, format="gztar", root_dir=TEMP_DIR, base_dir=".")

done()
