# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-08-28


# -------------- Required packages
import os
import h5py
import getpass
import re
import time


# -------------- Declare function(s)
def copy_group(source_group, dest_group):
    '''TODO
    '''
    # Loop over all groups in the source (group):
    for key in source_group:
        if isinstance(source_group[key], h5py.Group):
            # Create a new group:
            new_group = dest_group.create_group(key)

            # Copy attributes to the new group:
            for attr_name, attr_value in source_group[key].attrs.items():
                new_group.attrs[attr_name] = attr_value

            # This group might contain sub-groups, therefore, re-run function:
            copy_group(source_group[key], new_group)

        elif isinstance(source_group[key], h5py.Dataset):
            # Copy dataset to the group:
            dest_group.create_dataset(
                key,
                data=source_group[key][...],
                dtype=source_group[key].dtype,
                compression=source_group[key].compression,
                chunks=source_group[key].chunks
            )

            # Copy attributes to the group:
            for attr_name, attr_value in source_group[key].attrs.items():
                dest_group[key].attrs[attr_name] = attr_value

    return None


def copy_hdf5_file(source_file, destination_file):
    '''TODO
    '''
    with h5py.File(source_file, 'r') as src:
        with h5py.File(destination_file, 'w') as dest:
            # Copy attributes from the source file to the destination file:
            for attr_name, attr_value in src.attrs.items():
                dest.attrs[attr_name] = attr_value

            # Recursively copy all groups and datasets:
            copy_group(src, dest)

            # Time of snapshot:
            snap_time = src['Header'].attrs['Time']
            # dest['Header'].attrs['Time'] = 0

    return snap_time


def do_sweep_post_process(source_dir, snap_range):
    '''TODO

    Args:
        source_dir (str): TODO
    '''
    print('\nWelcome to this Arepo Sweep post-processing script ' +
          '(built using Vatpy)')
    print(f'  * Directory of snapshots to post-process: {source_dir}')
    print('  * Range of snapshots to post-process: '
          f'({snap_range[0]}, {snap_range[1]})')
    print('\nStarting the post-processing routine...\n')

    # Get current working directory:
    cwd = os.getcwd()

    # Check if neccessary files are present:
    print('INIT: Checking directory for requirred files needed to ' +
          'post-process snapshots')
    terminate = '(terminating this script)\n\nTerminating...\n\n'
    if not os.path.isfile('Arepo'):
        return print(f'ERROR: Missing Arepo file {terminate}')
    if not os.path.isfile('job.sh'):
        return print(f'ERROR: Missing job.sh file {terminate}')
    if not os.path.isfile('param.txt'):
        return print(f'ERROR: Missing param.txt file {terminate}')
    if not os.path.isfile('TreeCol_lookup.dat'):
        return print(f'ERROR: Missing TreeCol_lookup.dat file {terminate}')

    # Check if post-process directory and snapshots already exist(s):
    print('INIT: Checking if post-processed snapshots already exist in ' +
          'POSTPROCESSEDSNAPSHOTS')
    snap_start = int(snap_range[0])
    snap_end = int(snap_range[1])
    if os.path.isdir('POSTPROCESSEDSNAPSHOTS'):
        snap_str = '000'[:3-len(str(snap_start))] + str(snap_start)
        while os.path.isfile(f'POSTPROCESSEDSNAPSHOTS/snap_{snap_str}.hdf5'):
            if snap_start == snap_end:
                print('ENDING: All snapshots are already post-processed!\n')
                return 0
            else:
                snap_start += 1
                snap_str = '000'[:3-len(str(snap_start))] + str(snap_start)
        print(f'INIT: Found {snap_start - int(snap_range[0])} already ' +
              'post-processed snapshots')
    else:
        os.system('mkdir POSTPROCESSEDSNAPSHOTS')

    # Find user name:
    username = getpass.getuser()

    # Find job name:
    with open(cwd + '/job.sh', 'r') as file:
        for line in file:
            match_job = re.search(r'--job-name (\S+)', line)
            match_sweep = re.search(r'NUM_SWEEP_ITERATIONS=(\S+)', line)
            match_snap = re.search(r'NUM_ITERATIONS_PER_SNAPSHOT=(\S+)', line)
            if match_job:
                jobname = match_job.group(1)
            if match_sweep:
                num_sweep_iterations = int(match_sweep.group(1))
            if match_snap:
                num_iterations_per_snapshot = int(match_snap.group(1))

    print('INIT: Information to track when submitting jobs to SLURM:' +
          f'\n  USER: {username}\n  JOB-NAME: {jobname}')

    # Loop over remaining snapshots at source to post-process:
    for s in range(snap_start, snap_end+1):
        # Copy the source snapshot to the current working directory:
        source_file = (source_dir + '/snap_' + '000'[:3-len(str(s))] +
                       str(s) + '.hdf5')
        destination_file = cwd + '/postprocess.hdf5'
        print(f'PREP: Copying data from snapshot {s} to postprocess.hdf5')
        snap_time = copy_hdf5_file(source_file=source_file,
                                   destination_file=destination_file)

        # Read the default parameter file:
        with open('param.txt', 'r') as param:
            param_lines = param.readlines()

        for line in param_lines:
            if line.strip().startswith('MaxSizeTimestep'):
                delta_time = float(line.split()[-1])

        # Copy all parameters to a new parameter file,
        # but change TimeBegin and TimeMax:
        print('PREP: Creating a temporary parameter file')
        with open('param_postprocess.txt', 'w') as param_postprocess:
            for line in param_lines:
                if line.strip().startswith('MaxSizeTimestep'):
                    new_timestep_max = delta_time / num_sweep_iterations
                    param_postprocess.write(
                        f'MaxSizeTimestep    {new_timestep_max}\n')
                else:
                    param_postprocess.write(line)

        # Launch job:
        print(f'LAUNCH: About to post-process snapshot {s}')
        os.system('sbatch -Q job.sh')

        # Automatically check if the job is still running or not:
        time_start = time.time()
        while os.popen('squeue -u ' + username + ' -h -n ' + jobname +
                       ' -o "%.8j"').read():
            time.sleep(1)
            print(f'  ELAPSED TIME: {int(time.time() - time_start)} s       ',
                  end='\r')

        time_end = time.time()
        duration = int(time_end - time_start)

        # Check if the final output file exists:
        final = num_sweep_iterations // num_iterations_per_snapshot
        f_str = '000'[:3-len(str(final))] + str(final)
        if not os.path.isfile(f'OUTPUT/snap_{f_str}.hdf5'):
            print(f'WARNING: Snapshot {s} post-processed without generating ' +
                  f'a final output (took {duration} s)')
            print('ACTION: Please try to re-run the post-processing routine ' +
                  'using a different number of tasks')
            with open('SWEPPP-WARNINGS.txt', 'a') as warn:
                warn.write(f'Snapshot {s} post-processed without generating ' +
                           'a final output\n')
        else:
            print(f'COMPLETED: Snapshot {s} post-processed succesfully ' +
                  f'(took {duration} s)       ')

            # Move the final output to another directory:
            print('OUTPUT: Moving the final output to POSTPROCESSEDSNAPSHOTS')
            s_str = '000'[:3-len(str(s))] + str(s)
            os.system(f'mv OUTPUT/snap_{f_str}.hdf5 POSTPROCESSEDSNAPSHOTS/' +
                      f'snap_{s_str}.hdf5')

            # Update the time of the post-processed snapshot:
            file = f'POSTPROCESSEDSNAPSHOTS/snap_{s_str}.hdf5'
            with h5py.File(file, 'r+') as f:
                f['Header'].attrs['Time'] = snap_time

        print(f'DONE: Post-processing of snapshot {s} completed')

    print('ENDING: All snaphots post-processed!\n')

    return 0

# -------------- End of file
