#!/usr/bin/env bash

# adapted from: https://github.com/cbg-ethz/V-pipe

# defaults
PREFIX=$(pwd)
WORKDIR=

cat <<EOF

██╗    ██╗███████╗██╗    ██████╗ ██╗██████╗ ███████╗██╗     ██╗███╗   ██╗███████╗                        
██║    ██║██╔════╝██║    ██╔══██╗██║██╔══██╗██╔════╝██║     ██║████╗  ██║██╔════╝                        
██║ █╗ ██║█████╗  ██║    ██████╔╝██║██████╔╝█████╗  ██║     ██║██╔██╗ ██║█████╗                          
██║███╗██║██╔══╝  ██║    ██╔═══╝ ██║██╔═══╝ ██╔══╝  ██║     ██║██║╚██╗██║██╔══╝                          
╚███╔███╔╝██║     ██║    ██║     ██║██║     ███████╗███████╗██║██║ ╚████║███████╗                        
 ╚══╝╚══╝ ╚═╝     ╚═╝    ╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝╚═╝╚═╝  ╚═══╝╚══════╝                        
                                                                                                         
 █████╗ ██╗   ██╗████████╗ ██████╗ ██╗███╗   ██╗███████╗████████╗ █████╗ ██╗     ██╗     ███████╗██████╗ 
██╔══██╗██║   ██║╚══██╔══╝██╔═══██╗██║████╗  ██║██╔════╝╚══██╔══╝██╔══██╗██║     ██║     ██╔════╝██╔══██╗
███████║██║   ██║   ██║   ██║   ██║██║██╔██╗ ██║███████╗   ██║   ███████║██║     ██║     █████╗  ██████╔╝
██╔══██║██║   ██║   ██║   ██║   ██║██║██║╚██╗██║╚════██║   ██║   ██╔══██║██║     ██║     ██╔══╝  ██╔══██╗
██║  ██║╚██████╔╝   ██║   ╚██████╔╝██║██║ ╚████║███████║   ██║   ██║  ██║███████╗███████╗███████╗██║  ██║
╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝ ╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝╚═╝  ╚═╝
                                                                                             
EOF
sleep 1.0

# Helper
fail() {
	printf '\e[31;1mArgh: %s\e[0m\n'	"$1"	1>&2
	[[ -n "$2" ]] && echo "$2" 1>&2
	exit 1
}

oops() {
	printf '\e[33;1m %s\e[0m\n'	"$1"	 1>&2
	sleep 0.5
}

title() {
	printf '\e[34;1m======================\n%s\n======================\e[0m\n\n'	"$1"
	sleep 0.5
}

message() {
	printf '\e[32;1m%s\t%s\e[0m\n' "$1" "$2"
	sleep 0.5
}

status() {
	printf '\e[36;1m%s\e[0m\n'	"$1"
	sleep 0.5
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
-p PREFIX    prefix directory under which to install wfi
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
	status "No directory prefix specified, wfi will be installed in the current directory"
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
	if [[ -e $PROFILE ]] &&  grep -qH 'conda initialize\|CONDA\|module \(add\|load\)' $PROFILE; then
		ENVIRONMENTWARNING=1
		INSTALL_CONDA=0
	fi
done

if [[ -n "$CONDA_PREFIX" ]]; then
	message 'CONDA_PREFIX environment variable set'
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

	message 'installing Miniconda3'

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
	message 'Conda is already installed.'
	source $CONDA_PREFIX/etc/profile.d/conda.sh
	conda activate base
fi

status 'Checking Mamba installation'
if [[ $INSTALL_CONDA -eq "1" ]]; then
	if [[ -x $(which mamba) ]] && mamba -V 2>/dev/null >/dev/null; then
		message 'Wow mamba is installed! We got a hacker over here'
	else
		message 'Installing mamba, this might be quick or slow. Wait.'
		conda install --yes -n base -c conda-forge mamba
	fi
	conda config --add channels bioconda
	conda config --add channels conda-forge
fi

message 'Creating conda environment "wfi" and installing dependancies'
DEPEND=("r-ggplot2" "r-dplyr" "r-tidyr" "r-cowplot" "r-gridExtra" "r-optparse" "r-furrr")
ENVS=$(conda env list | grep 'wfi' )
if [[ -n $ENVS ]]; then
	conda activate wfi
	message "conda enviornment 'wfi' was detected."
	oops "Skipping conda and package installation; this assumes all dependancies are installed."
else 
	if mamba create -n wfi --yes -c agbiome -c anaconda python=3.6 r-forge snakemake-minimal cutadapt biopython bbtools $GIT ${DEPEND[@]} ; then
		message "Success! Conda environment created and deps installed."
	else
		fail "There was an error or conflict installing dependancies. Contact the administrator."
		exit 1
	fi
fi;


##############
#     wfi    #
##############

status 'Installing wfi'
message 'Attempting to download latest release of wfi'
LATEST_URL=$(curl -s https://api.github.com/repos/ammaraziz/wfi/releases | grep "browser_download_url" | cut -d '"' -f 4 | head -n 1)
NAME_WFI=$(basename $LATEST_URL)


wget -q $LATEST_URL -P $PREFIX
if [[ $? -ne 0 ]]; then
	oops "wget failed, trying --no-check-certificate"
	curl -s https://api.github.com/repos/ammaraziz/wfi/releases | grep "browser_download_url" | cut -d '"' -f 4 | head -n 1 | wget --no-check-certificate -q -	
fi
message 'Latest version downloaded, attempting to uncompress and setup wfi pipline.'

mkdir wfi_latest
tar -xf $PREFIX/$NAME_WFI -C wfi_latest
unzip -q wfi_latest/bin/flu-amd.zip -d wfi_latest/bin/
bash wfi_latest/tools/mod_init.sh -i wfi_latest/bin/flu-amd/IRMA_RES/modules/FLU/init.sh

message "wfi has finished installing!"

###############
#             #
#   Checks    #
#             #
###############

status 'Checking total Installation'
conda activate wfi

if [[ -x $(which conda) ]] && conda --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	CONDA=0
	message 'found conda!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'conda is missing, install manually.'
	CONDA=1
fi;

if [[ -x $(which snakemake) ]] && snakemake --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	SNAKE=0
	message 'found snakemake!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'snakemake is missing, install manually: conda install snakemake-minimal'
	SNAKE=1
fi;

if [[ -x $(which cutadapt) ]] && cutadapt --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	CUTADAPT=0
	message 'found cutadapt!'
else # out-of-the box Mac OS X and some Linux dockers
	oops 'cutadapt is missing, install manually: conda install {package}'
	CUTADAPT=1
fi;

if [[ -x $(which bbduk.sh) ]] && bbduk.sh --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	BBDUK=0
	message 'found bbduk!'
else
	oops 'bbduk is missing, install manually.'
	BBDUK=1
fi;

RPACKAGES=$(Rscript <(echo 'is_inst <- function(pkg) {nzchar(system.file(package = pkg))}
	p = c("ggplot2", "dplyr", "tidyr", "cowplot", "gridExtra", "furrr", "optparse")
	s = sum(unlist(lapply(p, is_inst)))
	if (s == 5) {write("0", stdout())} else {write("1", stdout())}'))

if [[ $RPACKAGES -eq "0" ]]; then
	message 'All R packages Installed!'
else
	oops 'Some or all R packages are missing, install manually via running R interactively'
fi;

# final check
if [[ $CONDA -eq "0" && $SNAKE -eq "0" && $CUTADAPT -eq "0" && $BBDUK -eq "0" && $RPACKAGES -eq "0" ]]; then
	title 'Installation of wfi complete! Remember to `conda activate wfi` before running the pipeline'
	exit 0
else
	oops 'Error during installation check. Look above'
	exit 1
fi