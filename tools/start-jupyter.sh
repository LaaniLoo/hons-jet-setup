#!/usr/bin/env bash

ncpus=4
notebook_dir='/scratch/${USER}/phd'
job_length=24

show_help()
{
    echo "-c = connect to an existing jupyter session"
    echo "-n = don't kill jupyter notebook when ssh tunnel closes"
    echo "-d = kill any existing jupyter session"
    echo "-l = a specific local port to use"
}

set_run_script()
{
    run_script=$(cat <<EOF
#!/bin/bash -l

#PBS -lselect=ncpus=${ncpus}:mpiprocs=${ncpus}
#PBS -lwalltime=${job_length}:00:00
#PBS -N jupyter
#PBS -j eo
#PBS -e notebook_log
conda activate analysis
module load intel impi hdf5 sz
port=\\\$(shuf -n 1 -i 9000-9300)
# notebook directory
notebook_dir=${notebook_dir}

# set the runtime directory
export JUPYTER_RUNTIME_DIR=~/.local/share/jupyter/runtime

echo \\\${port}
echo \\\${notebook_dir}
echo \\\$(hostname -s)

jupyter lab --ip=\\\$(hostname -s) --port=\\\${port} --notebook-dir=\\\${notebook_dir} --no-browser
EOF
)
}

check_for_jupyter()
{
    notebook_count=$(ssh kunanyi 'shopt -u nullglob; ls ~/.local/share/jupyter/runtime/nbserver*.json 2> /dev/null | wc -l')


    if [ "$notebook_count" -ne "0" ]; then
        if [ -n "$notebook_jobs" ]; then
            echo "The following jupyter job(s) are already running for your user:"
            echo "${notebook_jobs}"
            echo ""
            echo "Please either connect to an existing session with start-jupyter -c, or kill all existing sessions on kunanyi before starting a new one"
            echo "Running start-jupyter -d will kill existing jupyter jobs"
        else
            echo "${notebook_count} existing session files found, but no notebook jobs are running."
            echo "Try clearing all files in the ~/.local/share/jupyter/runtime directory on kunanyi"
        fi
        exit 1
    fi
}

start_jupyter() 
{
    echo "Submitting jupyterlab job on kunanyi"
    # notebook_job=$(ssh kunanyi 'qsub ~/jupyter-notebook.sh')
    notebook_job=$(ssh kunanyi bash -c "qsub - <<EOF
${run_script}
EOF")
    return
}

get_notebook_config()
{
    ssh kunanyi 'until compgen -G "$HOME/.local/share/jupyter/runtime/nbserver*.json"; do sleep 1; done' > /dev/null
    notebook_config=$(ssh kunanyi 'cat $(ls -rt ~/.local/share/jupyter/runtime/nbserver*.json | tail -n1)')
    notebook_port=$(echo $notebook_config | jq -r .port)
    notebook_hostname=$(echo $notebook_config | jq -r .hostname)
    notebook_url=$(echo $notebook_config | jq -r .url)
    notebook_token=$(echo $notebook_config | jq -r .token)
    # local_port=$(gnushuf -n 1 -i 9000-9400)
    return
}

get_notebook_jobs()
{
    notebook_jobs=$(ssh kunanyi 'qstat -wau "$USER" | grep "jupyter" | cut -d " " -f1')
    notebook_jobs_array=($notebook_jobs)
}

start_ssh_tunnel()
{
    echo "Connecting to jupyterlab on kunanyi"
    echo "Connection url is http://localhost:${local_port}/?token=${notebook_token}"
    ssh -L ${local_port}:${notebook_hostname}:${notebook_port} kunanyi -N
}

delete_notebook_job()
{

    echo "Deleting ${notebook_job} from kunanyi"
    ssh kunanyi qdel ${notebook_job}
}

delete_notebook_jobs()
{
    if [ -n "$notebook_jobs" ]; then
        # iterate over jobs
        for job_id in "${notebook_jobs_array[@]}"; do
            echo "Deleting ${job_id} from kunanyi"
            ssh kunanyi qdel ${job_id}
        done
    else
        echo "There are no jupyterlab instances running"
    fi
}

gnushuf()
{
    if command -v gshuf >/dev/null 2>&1; then
        gshuf "$@"
    else
        shuf "$@"
    fi
}

### MAIN ###

OPTIND=1

# no_exit=0
check_for_existing_session=1
start_jupyter_on_server=1
start_ssh_tunnel_on_client=1
stop_jupyter_on_server=1
local_port=$(gnushuf -n 1 -i 9000-9400)

while getopts "h?ncdl:" opt; do
    case "$opt" in 
        h|\?)
            show_help
            exit 0
            ;;
        l) local_port=${OPTARG}
            ;;
        n) stop_jupyter_on_server=0
            ;;
        c) 
            start_jupyter_on_server=0;
            start_ssh_tunnel_on_client=1;
            stop_jupyter_on_server=0;
            ;;
        d)
            start_jupyter_on_server=0;
            start_ssh_tunnel_on_client=0;
            stop_jupyter_on_server=1;
            ;;
    esac
done

shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift

if [ "$start_jupyter_on_server" -eq "1" ]; then
    if [ "$check_for_existing_session" -eq "1" ]; then
        get_notebook_jobs
        check_for_jupyter
    fi
    set_run_script
    start_jupyter
fi

if [ "$start_ssh_tunnel_on_client" -eq "1" ]; then
    if [ "$start_jupyter_on_server" -eq "1" ]; then
        echo "Waiting for job to start..."
    fi
    get_notebook_config;
    start_ssh_tunnel;

    echo "Closing connection"
fi

if [ "$stop_jupyter_on_server" -eq "1" ]; then
    if [ "$start_jupyter_on_server" -eq "1" ]; then
        delete_notebook_job
    else
        get_notebook_jobs
        delete_notebook_jobs
    fi
fi
