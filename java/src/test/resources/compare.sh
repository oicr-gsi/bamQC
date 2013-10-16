#!/bin/bash

usage() {
cat <<EOF
Usage for $0
$0 [-v] [-h] [-e {directory|file}] directory|file
-v  Verbose mode
-h  Show this usage menu
-e  Expected output
$verbose
EOF
exit 1
}

error() {
    echo "$2" 1>&2
    exit $1
}


# Get all arguments, 1st argument is not the script name
args=("$@")

# Options
verbose=false
expected=
input=

optspec=":vhe:"
while getopts "$optspec" OPT; do
	case "${OPT}" in
 	v)
		verbose=true
		;;
	h) 	
		usage
		;;	
	e)
		if [ -e ${OPTARG} ]; then
			expected=${OPTARG}
		else
			error 1 "Expected directory/file not accessible"
		fi
		;;
	:)
   		echo "Option '-${OPTARG}' missing argument"
   		;;	
	*)
		if [ "$OPTERR" != 1 ] || [ "${optspec:0:1}" = ":" ]; then
        	echo "Unknown option: '-${OPTARG}'" >&2
        fi
		;;
	esac
done

if [[ ! -z "${args[$OPTIND-1]}" && -e "${args[$OPTIND-1]}" ]]; then
	input=${args[$OPTIND-1]}
else
	error 1 "Input directory/file missing"
fi

calculate_output=$(
	##RUN CALCULATE SCRIPT
	cd "$input"
	#echo "error" 1>&2
	find . -type f -exec md5sum {} +
	#exit 3
	##
)

#Get exit status from calculate script
calculate_exit_status=$?

if (( $calculate_exit_status!=0 )); then
	error "$calculate_exit_status" "calculation step failed"
fi

if [ -z "$expected" ]; then
	echo "$calculate_output"
else
	##RUN COMPARE SCRIPT
	diff <(echo "$calculate_output") "$expected"
	##
	
	compare_exit_status=$?
	if (( $compare_exit_status!=0 )); then
		error "$compare_exit_status" "comparision step failed"
	fi
fi