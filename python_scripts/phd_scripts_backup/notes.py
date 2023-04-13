"""
                LET'S MAKE A NOTES APP
"""

import argparse
import datetime
import glob
import os
import pathlib


import sys, tempfile, os
from subprocess import call

PARSER = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Reweighting Scripts", epilog=" ")

PARSER.add_argument("-name", type=str, default='note',
                    help='Generic note (default: %(default)s)', required=False)
PARSER.add_argument("-dir", type=str, default='General',
                    help='Generic note (default: %(default)s)', required=False)
PARSER.add_argument("-ext", type=str, default='.md',
                    help='Generic note (default: %(default)s)', required=False)
PARSER.add_argument("-meeting", type=str, default=None,
                    help='new Meeting (default: %(default)s)', required=False)
PARSER.add_argument("-log", action="store_true",
                    help='Analyse SWISH data (default: %(default)s)')


ARGS = PARSER.parse_args()
EDITOR = os.environ.get('EDITOR','vim') #that easy!

def quicknote():
    initial_message = b"" # if you want to set up the file somehow
    with tempfile.NamedTemporaryFile(suffix=".md") as tf:
        tf.write(initial_message)
        tf.flush()
        call([EDITOR, tf.name])

        # do the parsing with `tf` using regular File operations.
        # for instance:
        tf.seek(0)
        edited_message = tf.read()
        new_content = edited_message.decode("utf-8")

    return new_content

def create_note(directory, name, ext, cont=True):
    path = pathlib.Path(directory+'/'+name+ext)
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

def suffix(d):
    return 'th' if 11<=d<=13 else {1:'st',2:'nd',3:'rd'}.get(d%10, 'th')

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
            f.write("\n\n# {} {} LOG FILE\n\n"
                    .format(today.strftime("%B").upper(),
                            today.strftime("%Y")))
    new_lines = quicknote()
    print('Adding to Log File: {}'.format(logfile))
    try:
        with path.open(mode='a') as f:
            f.write("\n## {t1}{s} {t2}\n".format(t1=today.strftime("%A %d"),
                                    s=suffix(today.day),
                                    t2=today.strftime("%B %Y")))
            f.write(new_lines)
    except:
        print('Saving Log failed')
        sys.exit()
    print('Successfully edited Log File: {}'.format(logfile))

def meeting(name):
    today = datetime.datetime.now()
    logfile = today.strftime("%Y-%m-%d_{}.md".format(name))
    path = pathlib.Path('Meetings/'+logfile)
    with path.open(mode='w+') as f:
        f.write("# {n} \n### {t1}{s} {t2} \n\n## Participants"
                .format(n=name.replace("_"," "),
                        t1=today.strftime("%A %d"),
                        s=suffix(today.day),
                        t2=today.strftime("%B %Y"),
                        ))
    create_note('Meetings',logfile,'')

if __name__ == '__main__':
    #create_note('./NEW_FOLDER/', 'TEST', '.txt')
    if ARGS.log:
        log()
    elif ARGS.meeting is not None:
        meeting(ARGS.meeting)
