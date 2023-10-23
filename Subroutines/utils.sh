function on_cluster() {
	if [[ "${CLUSTER-empty}" == "GREENE" ]]; then
		echo TRUE
	else
		echo FALSE
	fi
}

function which_software() {	
	if [[ "$(on_cluster)" == "TRUE" ]]; then
		echo SINGULARITY
	else
		echo DOCKER
	fi
}

function load_modules() {
	if [[ "$(on_cluster)" == "TRUE" ]]; then
		module load freesurfer/6.0.0
                module load matlab/2021a
	fi
}

function singularity_pull() {
	umask u=rwx,g=rx,o=rx
	singularity pull --force docker://${1}
}

function container_pull() {
	if [[ $CONTAINER_SOFTWARE == DOCKER ]]; then
		docker pull $@
	else		
		singularity_pull $@
	fi
}

function container_run_simple() {
	if [[ $CONTAINER_SOFTWARE == DOCKER ]]; then
		docker run $@
	else
		singularity run $@
	fi
}

function container_run() {
	if [[ $CONTAINER_SOFTWARE == DOCKER ]]; then
		docker run \
		--user "$(id -u):$(id -g)" \
		--rm \
		-v $1 \
		-v $2 \
		$3 \
		$5
	else
		singularity run \
		-B $1 \
		-B $2 \
		$4 \
		$5
	fi
}

function fsLicensePath() {
# FreeSurfer license path:
#      We first check whether FREESURFER_LICENSE is an environmnetal variable
if [ -z "${FREESURFER_LICENSE+set}" ]
then fsLicense=${FREESURFER_HOME}/license.txt
else fsLicense="$FREESURFER_LICENSE"
fi
[ -r "$fsLicense" ] || {
    echo "FreeSurfer license (${fsLicense}) not found!"
    echo "You can set a custom license path by storing it in the environment variable FREESURFER_LICENSE"
    exit 1
}
echo $fsLicense
}
