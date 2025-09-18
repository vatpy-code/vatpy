# Descritopn: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-15


# -------------- Required packages
import os


# -------------- FilmMaker
class FilmMaker:
    '''TODO
    '''
    def __init__(self, film, snap_range):
        # Variables:
        self.film = film
        self.snap_range = snap_range

        print('\nWelcome to Vatpy FilmMaker')
        print(f'  * Selected film: {self.film}')

    ##########################################################################
    ##########################################################################
    def search_for_frames(self):
        '''TODO
        '''
        snap_start = int(self.snap_range[0])
        snap_end = int(self.snap_range[1])

        # Check if vframes directory exists:
        if os.path.isdir('./vframes'):
            # Check if film directory exists:
            if os.path.isdir(f'./vframes/{self.film}'):
                # Check how many frames already have been generated:
                fnr = '000'[:3-len(str(snap_start))] + str(snap_start)
                while os.path.isfile(f'./vframes/{self.film}/' +
                                     f'{self.film}_{fnr}.png'):
                    if snap_start == snap_end:
                        print('  * All frames already generated!')
                        return 1
                    else:
                        snap_start += 1
                        fnr = '000'[:3-len(str(snap_start))] + str(snap_start)
                print('  * Found ' +
                      f'{int(snap_start - int(self.snap_range[0]))} ' +
                      'already generated frames!')
            else:
                os.system(f'mkdir ./vframes/{self.film}')
        else:
            os.system('mkdir ./vframes')

        return snap_start, snap_end

    ##########################################################################
    ##########################################################################
    def generate(self):
        '''TODO
        '''
        # Snapshot range:
        snap_start = int(self.snap_range[0])
        snap_end = int(self.snap_range[1])

        # Check already generated frames:
        print('  * Searching for already generated frames...')
        search = self.search_for_frames()
        if search != 1:
            snap_start = search[0]
            snap_end = search[1]

            # Generate plots according to selected film:
            print('  * Starting to generate frames...')
            if self.film == 'noctua':
                from .noctua import plot_noctua
                plot_noctua(snap_start, snap_end)
            else:
                print('  * No matching film!')
                print('  * Ending the script\n')

        # ffmpeg:
        print('  * Running ffmpeg to generate film...')
        self.ffmpeg(snap_start, snap_end)

        print('  * Done!\n')

        return None

    ##########################################################################
    ##########################################################################
    def ffmpeg(self, start, end):
        '''TODO
        '''
        fps = 15
        res_width = 1920
        res_height = 1080
        vframes = end - start + 1

        f = os.system(f'ffmpeg -r {fps} -f image2 -y' +
                      f' -s {res_width}x{res_height}' +
                      f' -start_number {start}' +
                      f' -i {os.getcwd()}/vframes/{self.film}/' +
                      f'{self.film}_%03d.png' +
                      f' -vframes {vframes}' +
                      ' -vcodec libx264' +
                      ' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white"' +
                      ' -crf 25' +
                      ' -pix_fmt yuv420p' +
                      f' {os.getcwd()}/{self.film}.mp4' +
                      ' > ffmpeg.out 2> ffmpeg.err')
        if f == 0:
            print('  * Film generated successfully!')
        else:
            print('  * Error: Something went wrong, please check the ffmpeg' +
                  ' output/error file')

        return None

# -------------- End of file
