'''
Description: TODO
Authour(s): Jonathan Petersson
Last updated: 2025-05-16
'''


# -------------- Required packages
import os
import h5py
import time


# -------------- Declare function(s)
def do_meshrelax(ic, runs, jobdir, modify, modify_args, filename, savepath,
                 wait='manual'):
    '''
    Description: TODO
    '''
    print('  * Now entering the mesh relaxation routine:')

    # Move IC file to mesh relax directory:
    os.system(f'mv {ic} {jobdir}')

    # Mesh relaxtion:
    for i in range(1, runs+1):
        print(f'    - Mesh relaxation: {i}/{runs}')

        # -------------- Submission -------------- #
        homedir = os.getcwd()  # your home directory
        os.chdir(jobdir)  # change directory

        # Submit job:
        os.system('sbatch -Q job.sh')

        # Wait until mesh relaxation is done:
        if wait == 'auto':
            # Automatic waiting:
            print('    - Please Note: You are running mesh relaxation on' +
                  ' \'auto\', i.e. the code will automatically\n' +
                  '      check the queue for \'relax\' & remain idle until' +
                  ' it no longer runs')
            start = time.time()
            while os.popen('squeue -u jpeterss -h -n relax -o "%.8j"').read():
                time.sleep(1)
                print(f'    - Elapsed time: {int(time.time() - start)} s    ',
                      end='\r')
        else:
            # Manual waiting:
            print('    - Please Note: You are running mesh relaxation on' +
                  ' \'manual\',\n' +
                  '      i.e. please keep track of the submitted job manually')
            jobdone = 'n'
            while jobdone == 'n':
                jobdone = input('    - Is the job completed? (y/n): ')

        os.chdir(homedir)  # go back to home directory

        # -------------- Submission -------------- #
        print('    - Mesh relaxation done!')

        # Check if all runs of mesh relaxation are complated:
        if i < runs:
            # If all runs are NOT completed, do some modification(s)
            # to the IC file:
            print('    - Doing some preperations for next mesh relaxation...')
            os.system('mv {jobdir}/OUTPUT/snap_001.hdf5' +
                      f' {jobdir}/{filename}.hdf5')
            if modify:
                print('    - Doing some modifications of the last snapshot...')
                modify(ic=f'{jobdir}/{filename}.hdf5', **modify_args)
            print('    - Preperations done!')

        elif i == runs:
            # If all runs are completed, move the IC file to the savepath:
            os.system(f'mv {jobdir}/OUTPUT/snap_001.hdf5' +
                      f' {savepath}/{filename}.hdf5')
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                f['Header'].attrs['Time'] = 0

    print('  * All iterations of mesh relaxation completed! Congrats!')
    print(f'  * New IC file: {savepath}/{filename}.hdf5')

    return 0


# -------------- End of file
