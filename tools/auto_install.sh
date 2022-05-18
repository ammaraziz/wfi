#!/usr/bin/env bash

# adapted from: https://github.com/cbg-ethz/V-pipe

# defaults
PREFIX=$(pwd)
WORKDIR=

# Helper
fail() {
	printf '\e[31;1mArgh: %s\e[0m\n'	"$1"	1>&2
	[[ -n "$2" ]] && echo "$2" 1>&2
	exit 1
}

oops() {
	printf '\e[33;1mS: %s\e[0m\n'	"$1"	 1>&2
}

title() {
	printf '\e[34;1m======================\n%s\n======================\e[0m\n\n'	"$1"
}

message() {
	printf '\e[32;1m%s\t%s\e[0m\n' "$1" "$2"
}

status() {
	printf '\e[36;1m%s\e[0m\n'	"$1"
}

check_directory() {
	if [[ -d "$1" ]]; then
		if (( FORCE )); then
			rm -rf "$1"
		else
			fail "${2:-Directory} ${1} already exists" 'check and use -f if needed'
		fi
	fi
}

usage() {
echo "usage: $0 [options]
options:
-f           force overwriting directories
-p PREFIX    prefix directory under which to install V-pipe
             [default: current directory]
-w WORKDIR   create and populate working directory
-m           only minimal working directory
-h           print this help message and exit"
}


# parameters
while getopts ":fp:b:r:w:mh" opt; do
	case ${opt} in
		f)
			FORCE=1
		;;
		p)
			PREFIX="${OPTARG}"
		;;
		w)
			WORKDIR="${OPTARG}"
		;;
		m)
			MINIMAL='-m'
		;;
		h)
			usage
			exit 0
		;;
		\?)
			fail "Invalid option: ${OPTARG}"
		;;
		:)
			fail "Invalid option: ${OPTARG} requires an argument"
		;;
	esac
done
shift $((OPTIND -1))

if [[ -n "${PREFIX}" ]]; then
	oops " No directory prefix specified, wfi will be installed in the current directory"
fi

###################
#                 #
#   Basic stuff   #
#                 #
###################

DOWNLOAD=
if [[ -x $(which wget) ]] || wget --version; then # most Linux distros (Gentoo)
	DOWNLOAD='wget'
elif [[ -x $(which curl) ]] || curl --version; then # a few Linux distros (CentOS, Archlinux) and Mac OS X
	DOWNLOAD='curl -O'
else
	fail 'Please install either wget or curl'
fi;

# HACK if not available locally, we will install git from conda later on (together with the other conda packages: snakemake, etc.)
GIT=
if [[ -x $(which git) ]] && git --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	GIT=
else # out-of-the box Mac OS X and some Linux dockers
	oops 'git is missing, I will download it from conda'
	GIT=git
fi;

# CHECK having environment modifiers (such as conda or modules) is a *horrendously bad idea* that can cause hard to understand errors much later
ENVIRONMENTWARNING=
for PROFILE in $HOME/.bash_profile $HOME/.bashrc $HOME/.profile; do
	if [[ -e $PROFILE ]] &&  grep -H 'conda initialize\|CONDA\|module \(add\|load\)' $PROFILE; then
		ENVIRONMENTWARNING=1
		INSTALL_CONDA=0
	fi
done

if [[ -n "$CONDA_PREFIX" ]]; then
	echo 'CONDA_PREFIX environment variable set'
	ENVIRONMENTWARNING=1
	INSTALL_CONDA=0
fi

##################
#                #
#   Miniconda3   #
#                #
##################

if ((INSTALL_CONDA)); then
	MINICONDA=
	MINICONDAPATH="${PREFIX}/miniconda3"

	title 'installing Miniconda3'

	# Check if directory is free
	check_directory "${MINICONDAPATH}" 'Miniconda3 installation path'

	# Check OS for OS-Spefic Miniconda3 installer
	if [[ "$OSTYPE" == linux* ]]; then
		MINICONDA=Miniconda3-latest-Linux-x86_64.sh
	elif [[ "$OSTYPE" == darwin* ]]; then
		MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh
	else
		fail 'I cannot detect OS. Only Linux and Mac OS X are supported' 'manually override OSTYPE environment variable if needed'
	fi

	message 'Using installer:' ${MINICONDA}

	# Get and install Miniconda3
	mkdir -p ${PREFIX}
	cd ${PREFIX}
	[[ -f ${MINICONDA} ]] && rm ${MINICONDA}
	${DOWNLOAD} https://repo.anaconda.com/miniconda/${MINICONDA}
	sh ${MINICONDA} -b -p miniconda3
	# -b for batch (no question asked)

	. miniconda3/bin/activate
else
	title 'Conda is already installed!'
fi

# channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

title 'Creating conda environment "wfi" and installing dependancies'

DEPEND=(r-ggplot2 r-dplyr r-tidyr r-cowplot r-gridExtra r-optparse)

conda install --yes -n base -c conda-forge mamba

ENVS=$(conda env list | grep 'wfi' )
if [[ -n $ENVS ]]; then
	source activate wfi
	message "conda enviornment 'wfi' was detected. Skipping conda and package installation; this assumes all dependancies are installed! "
else 
	if mamba create -n wfi --yes python=3.6 r-forge snakemake-minimal $GIT $DEPEND cutadapt biopython; then
		message "Success! Conda environment created and deps installed."
	else:
		message "There was an error or conflict installing dependancies. Contact the administrator."
	fi
fi;


##############
#            #
#   wfi      #
#            #
##############

title 'Installing wfi'

message 'Attempting to download latest release of wfi'
# fix the below to /releases/latest when ready
curl -s https://api.github.com/repos/ammaraziz/wfi/releases | grep "browser_download_url" | cut -d '"' -f 4 | head -n 1 | wget -qi -
if [[ $? -ne 0 ]]; then
	oops "wget failed, trying --no-check-certificate"
	curl -s https://api.github.com/repos/ammaraziz/wfi/releases | grep "browser_download_url" | cut -d '"' -f 4 | head -n 1 | wget --no-check-certificate -qi -	
fi

# uncompress directory and install modules
tar -xvf *.tar.gz 
WFI_DIR=$(ls -d wfi*/)
unzip -q $WFI_DIR/bin/flu-amd.zip -d $WFI_DIR/bin/
cp -r $WFI_DIR/bin/custom_modules/RSV/ $WFI_DIR/bin/flu-amd/IRMA_RES/modules/RSV/
bash $WFI_DIR/tools/mod_init.sh -i $WFI_DIR/bin/flu-amd/IRMA_RES/modules/FLU/init.sh

messsage "wfi has finished installing"



echo $'\n'

###############
#             #
#   Checks    #
#             #
###############

title 'Checking wfi Installation'

if [[ -x $(which conda) ]] && conda --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	CONDA=1
	message 'found conda!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'conda is missing, install manually.'
	CONDA=0
fi;

if [[ -x $(which snakemake) ]] && snakemake --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	SNAKE=1
	message 'found snakemake!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'cutadapt is missing, install manually: conda install snakemake-minimal'
	SNAKE=0
fi;

if [[ -x $(which cutadapt) ]] && cutadapt --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	CUTADAPT=1
	message 'found cutadapt!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'cutadapt is missing, install manually: conda install {package}'
	CUTADAPT=0
fi;

RPACKAGES=$(Rscript <(echo 'is_inst <- function(pkg) {nzchar(system.file(package = pkg))}
	p = c("ggplot2", "dplyr", "stringr", "tidyr", "cowplot", "gridExtra")
	s = sum(unlist(lapply(p, is_inst)))
	if (s <= 5) {write("0", stdout())} else {write("1", stdout())}'))

if (( RPACKAGES )); then
	message 'All R packages Installed!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'Some or all R packages are missing, install manually: conda install {r-package}'
	exit
fi;

echo $'\n'
title 'Installation of wf completed'
exit 0
