import os

sequence_folder = '/home/adrian/projects/brca2/data/raw_fastq/haches/'

folders = os.listdir(sequence_folder)
folders.sort()
print(folders)

for folder in folders:
    path = sequence_folder + folder

    print('working with {}'.format(path))
    os.chdir(path)

    files = os.listdir(path)
    for file in files:
        os.system('zcat {} | wc -l'.format(file))
    print()
