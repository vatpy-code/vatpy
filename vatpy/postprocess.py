'''
Description: TODO
Last updated: 2025-08-28
Authour(s): Jonathan Petersson
'''


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
            # Copy attributes to the group:
            for attr_name, attr_value in source_group[key].attrs.items():
                dest_group[key].attrs[attr_name] = attr_value

            # Copy dataset to the group:
            dest_group.create_dataset(
                key,
                data=source_group[key][...],
                dtype=source_group[key].dtype,
                compression=source_group[key].compression,
                chunks=source_group[key].chunks
            )

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

    return None


def do_sweep_post_process(source_dir, snap_range):
    '''TODO

    Args:
        source_dir (str): TODO
    '''
    print('\nWelcome to this Arepo Sweep post-processing script ' +
          '(built using Vatpy)')
    print(f'INFORMATION: Directory of snapshots to post-process: {source_dir}')
    print('INFORMATION: Range of snapshots to post-process: '
          f'({snap_range[0]}, {snap_range[1]})')

    # Get current working directory:
    cwd = os.getcwd()

    # Check if neccessary files are present:
    print('INITIALISING: Checking directory for requirred files needed to ' +
          'post-process Arepo snapshots')
    if not os.path.isfile(f'{cwd}/Arepo'):
        return 1
    if os.path.isfile(f'{cwd}/job.sh'):
        return 1
    if os.path.isfile(f'{cwd}/param.txt'):
        return 1
    if os.path.isfile(f'{cwd}/TreeCol_lookup.dat'):
        return 1

    # Check if post-process directory and snapshots already exist(s):
    print('INITIALISING: Checking if post-processed snapshots already exist')
    snap_start = snap_range[0]
    snap_end = snap_range[1]
    if os.path.isdir(f'{cwd}/POSTPROCESSEDSNAPHOTS'):
        snap_str = '000'[:3-len(str(snap_start))] + str(snap_start)
        while os.path.isfile(f'{cwd}/POSTPROCESSEDSNAPHOTS/snap_{snap_str}'):
            snap_start += 1
            snap_str = '000'[:3-len(str(snap_start))] + str(snap_start)
        print(f'INITIALISING: Found {snap_start - snap_range[0]} already ' +
              'post-processed snapshots')
    else:
        os.makedir(f'{cwd}/POSTPROCESSEDSNAPHOTS')

    if snap_start == snap_end:
        print('  ENDING: All snapshots are already post-processed!')
        return 0

    # Find user name:
    username = getpass.getuser()

    # Find job name:
    with open(cwd + '/job.sh', 'r') as file:
        for line in file:
            jobmatch = re.search(r'--job-name (\S+)', line)
            if jobmatch:
                jobname = jobmatch.group(1)

    print('INITIALISING: About to post-process snapshots for user={username}' +
          f'using job-name={jobname}')

    # Loop over remaining snapshots at source to post-process:
    for s in range(snap_start, snap_end):
        # Copy the source snapshot to the current working directory:
        source_file = (source_dir + '/snap' + '000'[:3-len(str(snap_start))] +
                       str(snap_start) + '.hdf5')
        destination_file = cwd + '/postprocess.hdf5'
        copy_hdf5_file(source_file=source_file,
                       destination_file=destination_file)

        # Launch job:
        print('LAUNCH: submitting post-processing job for snapshot {s}')
        os.system('sbatch -Q job.sh')

        # 
        time_start = time.time()
        while os.popen('squeue -u ' + username + ' -h -n ' + jobname +
                       '-o "%.8j"').read():
            time.sleep(2)
            print(f'Elapsed time: {int(time.time() - time_start)} s    ',
                  end='\r')

    return 0

# -------------- End of file
