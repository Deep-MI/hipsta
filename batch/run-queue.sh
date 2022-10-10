#!/bin/bash

if [ $# -ne 8 ] ; then

    echo
    echo "-------------------"
    echo "run-queue.sh"
    echo "-------------------"
    echo
    echo "This is a script to run multiple processes in parallel in a multi-core environment."
    echo
    echo "The script will read a list of IDs and a list commands, which must "
    echo "correspond to each other. From a given starting position in these lists, "
    echo "the script will extract a given number of cases, split these into a given "
    echo "number of processing queues with one or more cases each, and execute these "
    echo "queues in parallel. Any terminal outputs of the command will be saved to a "
    echo "case-specific log-file in the subjects directory, for which a suffix can be "
    echo "supplied. The script will also generate output about which case is currently "
    echo "proessed in which queue, including their process IDs."
    echo
    echo "The script has been written for and tested with the Bourne-Again Shell, "
    echo "i.e. the 'bash' environment, and should be run from such an environment. "
    echo
    echo
    echo "usage: run-queue.sh <NCASES> <NQUEUES> <START> <SUBJECTS_DIR> <LIST_FILE> <CMD_FILE> <SUFFIX> <TIMEOUT>"
    echo
    echo "  NCASES          ... number of cases to run"
    echo "  NCQUEUES        ... number of queues to run in parallel"
    echo "  START           ... starting index (zero-based)"
    echo "  SUBJECTS_DIR    ... subjects directory"
    echo "  LIST_FILE       ... file with list of subject IDs"
    echo "  CMD_FILE        ... file with list of commands to execute, one per subject ID"
    echo "                      must have same order and same number of lines as the list"
    echo "                      of subject IDs "
    echo "  SUFFIX          ... suffix appended to logfiles, e.g. to distinguish different"
    echo "                      analysis runs"
    echo "  TIMEOUT         ... timeout in seconds after which the processing of the current"
    echo "                      case is terminated (timeout 0 means no timeout)"
    echo
    echo "The recommended way to run this script is to use the 'nohup' command to detach "
    echo "it from a particular terminal, and to let it run in the background using the '&'"
    echo "character:"
    echo
    echo "  nohup run-queue.sh <...> > logfile.log &"
    echo
    echo "Note that it may make sense to limit the number of cores per process to one when"
    echo "using this program, e.g. by setting environment variables like so:"
    echo
    echo "  export OMP_NUM_THREADS=1"
    echo "  export ITK_NUM_THREADS=1"
    echo

    exit 0

fi

# ------------------------------------------------------------------------------
# subfunction
# ------------------------------------------------------------------------------

run_queue_indices() {

    k=0

    for i in ${CMDQUEUE[@]} ; do

        {

            sleep 2 ;

            date > ${SUBJECTS_DIR}/${LISTARRAY[$((i))]}${SFX}.log ;

            echo "Queue ${IDXQUEUE}, item ${k}: ${LISTARRAY[$((i))]} : started at "`date`

            # don't forget that we are using zero-based indicis originally
            if [ $TIMEOUT -gt 0 ] ; then
                timeout $TIMEOUT ` cat $CMD | head -n $((i+1)) | tail -n 1 ` >> ${SUBJECTS_DIR}/${LISTARRAY[$((i))]}${SFX}.log 2>&1 ;
            else
                ` cat $CMD | head -n $((i+1)) | tail -n 1 ` >> ${SUBJECTS_DIR}/${LISTARRAY[$((i))]}${SFX}.log 2>&1 ;
            fi

            if [ $? -eq 0 ] ; then
                echo "Queue ${IDXQUEUE}, item ${k}: ${LISTARRAY[$((i))]} : finished at "`date` ;
            else
                echo "Queue ${IDXQUEUE}, item ${k}: ${LISTARRAY[$((i))]} : exited with code $? at "`date` ;
            fi

            date >> ${SUBJECTS_DIR}/${LISTARRAY[$((i))]}${SFX}.log ;

        }

        k=$((k+1))

    done

}


# ------------------------------------------------------------------------------
# main script
# ------------------------------------------------------------------------------

# get input arguments

NCASES=$1
NQUEUES=$2
START=$3
SUBJECTS_DIR=$4
LIST=$5
CMD=$6
SFX=$7
TIMEOUT=$8

# set num threads

export OMP_NUM_THREADS=1
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

# calculate some numbers

NLIST=`cat $LIST | wc -w`

NCASESPERQUEUE=`echo "scale=0 ; ${NCASES}/${NQUEUES} + ((${NCASES}%${NQUEUES})>0)" | bc`

# check if dir exists

if [ ! -d $SUBJECTS_DIR ] ; then
    echo ""
    echo "Could not find subjects directory ${SUBJECTS_DIR}, exiting."
    echo ""
    exit 1
fi

# say something

echo
echo "Running ${NCASES} cases from ${LIST}, starting from ${START} and split into ${NQUEUES} queues."
echo "Starting at "`date`" with PID $$"
echo

# create array from list file

LISTARRAY=(`cat $LIST`)

# this is the main loop

IDXQUEUE=0 ; while [ $IDXQUEUE -lt $NQUEUES ] ; do

    if [ $((NCASES-IDXQUEUE*NCASESPERQUEUE)) -lt $NCASESPERQUEUE ] ; then

        NCASESPERCURRENTQUEUE=$((NCASES-IDXQUEUE*NCASESPERQUEUE))

    else

        NCASESPERCURRENTQUEUE=$NCASESPERQUEUE

    fi

    QUEUE=(`echo ${LISTARRAY[@]:$((START+IDXQUEUE*NCASESPERQUEUE)):${NCASESPERCURRENTQUEUE}}`)

    CMDQUEUE=(`seq -s ' ' $((START+IDXQUEUE*NCASESPERQUEUE)) $((START+IDXQUEUE*NCASESPERQUEUE+NCASESPERCURRENTQUEUE-1))`)

    run_queue_indices &

    echo "Queue $IDXQUEUE (PID $!) : ${QUEUE[@]}"

    IDXQUEUE=$((IDXQUEUE+1))

done

echo

wait

echo
echo "Finished at "`date`

#
