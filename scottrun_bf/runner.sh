#!/bin/bash
case $# in
        0|1)
                echo "Usage: runner.sh ifirst ifinal  // runs from i=ifirst to <=ifinal";
        exit 1 ;;
        2)
                firsti=$1
                lasti=$2
                for ((ii=${firsti};ii<=${lasti};ii++))
                do
                        echo "____________ msubal for run number " ${ii} ______________;
                        msubal pars${ii};
                done
esac