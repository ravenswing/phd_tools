"----------    VIMRC    ----------

" rebind the map leader to ,
let mapleader=","

"----------    PLUGINS    ----------
call plug#begin('~/.vim/plugged')

" nord colour scheme
Plug 'arcticicestudio/nord-vim'
let g:nord_italic = 1
let g:nord_underline = 1
let g:nord_italic_comments = 1

" vimtex for LaTeX support and compiling
Plug 'lervag/vimtex'
let g:tex_flavor='latex'
let g:vimtex_view_method='zathura'
let g:vimtex_quickfix_mode=0
set conceallevel=1
let g:tex_conceal='abdmg'

" NERDTree for file navigation
Plug 'scrooloose/nerdtree'

" Syntastic for syntax checking on load and save
Plug 'vim-syntastic/syntastic'
let g:syntastic_always_populate_loc_list = 1
let g:syntastic_check_on_open = 1
let g:syntastic_check_on_wq = 0
let g:syntastic_error_symbol = "✗"

" PEP8 Python format checking
Plug 'nvie/vim-flake8'

call plug#end()

" enable colour scheme
if has("termguicolors")
    set termguicolors
endif

colorscheme nord

set tabstop=4           " a tab is four spaces
set softtabstop=4
set shiftwidth=4        " number of spaces to use for autoindenting
set shiftround          " use multiple of shiftwidth when indenting with '<' and '>'
set expandtab

set colorcolumn=80      " add line to show 80 character 'limit'
set autoindent          " always set autoindenting on
set copyindent          " copy the previous indentation on autoindenting

set number relativenumber
set ruler
set showcmd 
filetype indent on
set wildmenu 
set showmatch   

set backspace=eol,start,indent " allow backspacing over everything in insert mode

set incsearch
set hlsearch 
set ignorecase          " ignore case if search pattern is all lowercase
set smartcase           " case-sensitive otherwise

set history=1000        " remember more commands and search history
set undolevels=1000     " use many muchos levels of undo
set title               " change the terminal's title
set visualbell          " don't beep
set noerrorbells        " don't beep

set undofile            " adds an undo file for unlimited undos
set gdefault            " set find/replace to global PER LINE as default

syntax enable 
let python_highlight_all=1

" remap to reduce keystrokes for all commands
nnoremap ; :
" move vertically by visual line
nnoremap j gj
nnoremap k gk

" move to beginning/end of line
nnoremap B ^
nnoremap E $

" $/^ doesn't do anything
nnoremap $ <nop>
nnoremap ^ <nop>

"split navigations
nnoremap <C-h> <C-w><C-h>
nnoremap <C-j> <C-w><C-j>
nnoremap <C-k> <C-w><C-k>
nnoremap <C-l> <C-w><C-l>

map <C-n> :NERDTreeToggle<CR>

" Use Q for formatting the current paragraph (or selection)
vmap Q gq
nmap Q gqap

" Clear search highlighting with ,/
nmap <silent> <leader>/ :nohlsearch<CR>

" Save all files when tabbing away - save on losing focus
au FocusLost * :wa

" Set F2 to open paste mode
nnoremap <F2> :set invpaste paste?<CR>
set pastetoggle=<F2>
set showmode

" Add spell-check and auto-correct with Ctrl+L
autocmd BufRead,BufNewFile *.md setlocal spell
autocmd BufRead,BufNewFile *.txt setlocal spell
autocmd BufRead,BufNewFile *.tex setlocal spell
set spelllang=nl,en_gb
inoremap <C-l> <c-g>u<Esc>[s1z=`]a<c-g>u

silent autocmd BufNewFile,BufRead *.py
"    \ set tabstop=4
"    \ set softtabstop=4
"    \ set shiftwidth=4
"    \ set textwidth=79
"    \ set expandtab
    \ set fileformat=unix
    \ set encoding=utf-8

"au BufRead,BufNewFile *.py,*.pyw,*.c,*.h match BadWhitespace /\s\+$/

" Always show statusline
"set laststatus=2

