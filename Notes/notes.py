"""
                LET'S MAKE A NOTES APP
"""

import argparse
import datetime
import glob
import os
import pathlib



PARSER = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Reweighting Scripts", epilog=" ")

PARSER.add_argument("-name", type=str, default='note',
                    help='Generic note (default: %(default)s)', required=False)
PARSER.add_argument("-dir", type=str, default='General',
                    help='Generic note (default: %(default)s)', required=False)
PARSER.add_argument("-ext", type=str, default='.md',
                    help='Generic note (default: %(default)s)', required=False)
PARSER.add_argument("-log", action="store_true",
                    help='Analyse SWISH data (default: %(default)s)')


ARGS = PARSER.parse_args()

def create_note(directory, name, ext, cont=True):
    path = pathlib.Path(directory+name+ext)
    if not path.parent.exists():
        path.parent.mkdir(parents=True)
    command = "gnome-terminal --command='vim {}'".format(path)
    if cont:
        command += ' &'
    print(command)
    os.system(command)

def unique_path(directory, name_pattern):
    counter = 0
    while True:
        counter += 1
        path = directory / name_pattern.format(counter)
        if not path.exists():
            return path
    # Usage:
    #path = unique_path(pathlib.Path.cwd(), 'test{:03d}.txt')

def unnamed_note(directory=None):
    if directory is not None:
        path = unique_path(pathlib.Path(directory), 'Note_{:02d}.md')
    else:
        path = unique_path(pathlib.Path.cwd(), 'Note_{:02d}.md')
    command = "gnome-terminal --command='vim {}' &".format(path)
    os.system(command)

def log():
    today = datetime.datetime.now()
    logfile = today.strftime("%Y_%B.md")
    path = pathlib.Path('Research_Logs/'+logfile)
    if not path.exists():
        with path.open(mode='w+') as f:
            f.write("\n\n### {} {} LOG FILE\n\n"
                    .format(today.strftime("%B").upper(),
                            today.strftime("%Y")))
    create_note('Research_Logs','log','.tmp', cont=False)
    print('THIS COMES AFTER EDITING')
    print(logfile)


if __name__ == '__main__':
    #create_note('./NEW_FOLDER/', 'TEST', '.txt')
    log()
