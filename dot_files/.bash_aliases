#BASH_ALIASES

# PhD specific
alias mountwd='sudo mount /dev/sdb1 /workdrive/'
alias gpfs='sudo mount -t cifs -o domain=ad.ucl.ac.uk,user=zcqsrev,uid=rhys //live.rd.ucl.ac.uk/ritd-ag-project-rdprod-flger36 /home/rhys/Desktop/gpfs'
alias unmount='sudo umount /dev/sdb1'
alias wd='cd /media/rhys/ExtHD/Project/'
alias uclvpn='sudo openconnect vpn.ucl.ac.uk --csd-wrapper ~/.config/csd-wrapper.sh'
alias ib='ssh rhys@ib-server.chem.ucl.ac.uk'
alias mn='ssh cnio96742@mn2.bsc.es'
alias th='ssh zcqsrev@thomas.rc.ucl.ac.uk'
alias gr='ssh zcqsrev@grace.rc.ucl.ac.uk'
alias my='ssh zcqsrev@myriad.rc.ucl.ac.uk'
alias jd='ssh -l rxe15-ffg01 jade.hartree.stfc.ac.uk'
alias jd2='ssh -l cxe41-ffg01 jade.hartree.stfc.ac.uk'
alias ac='ssh -l rhys login.archer.ac.uk'
alias ac2='ssh -l rhyse login.archer.ac.uk'
alias amran='ssh -X amran@flg230.chem.ucl.ac.uk'
alias carol='ssh -X arc@flg243.chem.ucl.ac.uk' # use if out of office. 
alias pc='ssh -X rhys@flg234.chem.ucl.ac.uk' # use if out of office. 
alias el='ssh revans@ela.cscs.ch'

alias lab='wd && cd Notebooks/ && jupyter lab'
alias manual='zathura ~/Dropbox/coding_manual/Coding_Manual.pdf &'

# Programs
alias vi='vim'
alias pymol='cd ~/software/pymol; ./pymol'
alias spotify='spotify --force-device-scale-factor=2.0 &'
alias gmx='/usr/local/gromacs/bin/gmx'

# Navigation
alias ls='ls -F --color=auto --group-directories-first'
alias ll='ls -ltrh'
alias la='ls -alh'
alias l='ls'
alias lt='ls -ltrh'
alias ..="cd .."
alias fa="find . -name '*'"
alias histg="history | grep "

alias bashrc='vi -p ~/.bash_aliases ~/.bashrc'
alias resrc='source ~/.bashrc'

# Functions
mcd () {
    mkdir -p $1
    cd $1
}
cs () {
    cd $1
    l
}
backup () {
    cp -r .config/ ~/Dropbox/Backup/Code/DotFiles/
    cp -r .vim/ ~/Dropbox/Backup/Code/DotFiles/
    cp {.bashrc,.bash_aliases,.bash_profile,.vimrc} ~/Dropbox/Backup/Code/DotFiles/
    wd
    cp *.py ~/Dropbox/Backup/Code/
    cp *.sh ~/Dropbox/Backup/Code/
    cp -r Gromacs/* ~/Dropbox/Backup/Code/Gromacs/
    cp -r Scripts/* ~/Dropbox/Backup/Code/Scripts/
    cp -r Notebooks/* ~/Dropbox/Backup/Notebooks/
}

function extract {
 if [ -z "$1" ]; then
    # display usage if no parameters given
    echo "Usage: extract ."
 else
if [ -f $1 ] ; then
        # NAME=${1%.*}
        # mkdir $NAME && cd $NAME
        case $1 in
          *.tar.bz2) tar xvjf ../$1 ;;
          *.tar.gz) tar xvzf ../$1 ;;
          *.tar.xz) tar xvJf ../$1 ;;
          *.lzma) unlzma ../$1 ;;
          *.bz2) bunzip2 ../$1 ;;
          *.rar) unrar x -ad ../$1 ;;
          *.gz) gunzip ../$1 ;;
          *.tar) tar xvf ../$1 ;;
          *.tbz2) tar xvjf ../$1 ;;
          *.tgz) tar xvzf ../$1 ;;
          *.zip) unzip ../$1 ;;
          *.Z) uncompress ../$1 ;;
          *.7z) 7z x ../$1 ;;
          *.xz) unxz ../$1 ;;
          *.exe) cabextract ../$1 ;;
          *) echo "extract: '$1' - unknown archive method" ;;
        esac
else
echo "$1 - file does not exist"
    fi
fi
}

# Make some possibly destructive commands more interactive.
#alias rm='rm -i'
alias mv='mv -i'
alias mkdir="mkdir -pv"

# Make grep more user friendly by highlighting matches
# and exclude grepping through .svn folders.
alias grep='grep --color=auto --exclude-dir=\.svn'

#if [ -f /usr/local/lib/python2.7/dist-packages/powerline/bindings/bash/powerline.sh ]; then
 #   source /usr/local/lib/python2.7/dist-packages/powerline/bindings/bash/powerline.sh
#fi
