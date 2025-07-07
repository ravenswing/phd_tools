---------------------------- NVCHAD options --------------------------
vim.o.laststatus = 3
vim.o.showmode = false

vim.o.clipboard = "unnamedplus"
vim.o.cursorline = true
vim.o.cursorlineopt = "number"

vim.opt.fillchars = { eob = " " }
vim.o.ignorecase = true
vim.o.smartcase = true

-- Numbers
vim.o.numberwidth = 2

-- disable nvim intro
vim.opt.shortmess:append "sI"

vim.o.signcolumn = "yes"
vim.o.splitbelow = true
vim.o.splitright = true
-- vim.o.timeoutlen = 400
vim.o.undofile = true

-- go to previous/next line with h,l,left arrow and right arrow
-- when cursor reaches end/beginning of line
vim.opt.whichwrap:append "<>[]hl"

-- disable some default providers
vim.g.loaded_node_provider = 0
vim.g.loaded_python3_provider = 0
vim.g.loaded_perl_provider = 0
vim.g.loaded_ruby_provider = 0

-- add binaries installed by mason.nvim to path
local is_windows = vim.fn.has "win32" ~= 0
local sep = is_windows and "\\" or "/"
local delim = is_windows and ";" or ":"
vim.env.PATH = table.concat({ vim.fn.stdpath "data", "mason", "bin" }, sep) .. delim .. vim.env.PATH

--------------------- MINE ---------------------------

vim.opt.number = true
vim.opt.relativenumber = true
vim.opt.ruler = false

vim.opt.tabstop = 4
vim.opt.softtabstop = 4
vim.opt.shiftwidth = 4
vim.opt.expandtab = true

vim.opt.smartindent = true

vim.opt.wrap = false

vim.opt.swapfile = false
vim.opt.backup = false
vim.opt.undodir = os.getenv "HOME" .. "/.vim/undodir"
vim.opt.undofile = true

vim.opt.hlsearch = false
vim.opt.incsearch = true

vim.opt.termguicolors = true

vim.opt.scrolloff = 8
vim.opt.signcolumn = "yes"
vim.opt.isfname:append "@-@"

vim.opt.updatetime = 50

vim.opt.colorcolumn = "80"
