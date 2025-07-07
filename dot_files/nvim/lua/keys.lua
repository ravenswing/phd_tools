------------- NV Chad bindings --------------
vim.keymap.set("i", "<C-b>", "<ESC>^i", { desc = "move beginning of line" })
vim.keymap.set("i", "<C-e>", "<End>", { desc = "move end of line" })
vim.keymap.set("n", "<Esc>", "<cmd>noh<CR>", { desc = "general clear highlights" })

vim.keymap.set("n", "<leader>ch", "<cmd>NvCheatsheet<CR>", { desc = "toggle nvcheatsheet" })

vim.keymap.set({ "n", "x" }, "<leader>fm", function()
  require("conform").format { lsp_fallback = true }
end, { desc = "general format file" })

-- global lsp mappings
vim.keymap.set("n", "<leader>ds", vim.diagnostic.setloclist, { desc = "LSP diagnostic loclist" })

-- tabufline
-- vim.keymap.set("n", "<leader>b", "<cmd>enew<CR>", { desc = "buffer new" })
-- vim.keymap.set("n", "<tab>", function()
-- require("nvchad.tabufline").next()
-- end, { desc = "buffer goto next" })
-- vim.keymap.set("n", "<S-tab>", function()
-- require("nvchad.tabufline").prev()
-- end, { desc = "buffer goto prev" })
-- vim.keymap.set("n", "<leader>x", function()
-- require("nvchad.tabufline").close_buffer()
-- end, { desc = "buffer close" })

-- Comment
vim.keymap.set("n", "<leader>/", "gcc", { desc = "toggle comment" })
vim.keymap.set("v", "<leader>/", "gc", { desc = "toggle comment" })

-- nvimtree
vim.keymap.set("n", "<C-n>", "<cmd>NvimTreeToggle<CR>", { desc = "nvimtree toggle window" })
-- vim.keymap.set("n", "<leader>e", "<cmd>NvimTreeFocus<CR>", { desc = "nvimtree focus window" })

-- vim.keymap.set("n", "<leader>th", function()
-- require("nvchad.themes").open()
-- end, { desc = "telescope nvchad themes" })

-- terminal
vim.keymap.set("t", "<C-x>", "<C-\\><C-N>", { desc = "terminal escape terminal mode" })

-- new terminals
vim.keymap.set("n", "<leader>h", function()
  require("nvchad.term").new { pos = "sp" }
end, { desc = "terminal new horizontal term" })

vim.keymap.set("n", "<leader>v", function()
  require("nvchad.term").new { pos = "vsp" }
end, { desc = "terminal new vertical term" })

-- toggleable
vim.keymap.set({ "n", "t" }, "<A-v>", function()
  require("nvchad.term").toggle { pos = "vsp", id = "vtoggleTerm" }
end, { desc = "terminal toggleable vertical term" })

vim.keymap.set({ "n", "t" }, "<A-h>", function()
  require("nvchad.term").toggle { pos = "sp", id = "htoggleTerm" }
end, { desc = "terminal toggleable horizontal term" })

vim.keymap.set({ "n", "t" }, "<A-i>", function()
  require("nvchad.term").toggle { pos = "float", id = "floatTerm" }
end, { desc = "terminal toggle floating term" })

-- LSP
vim.keymap.set("n", "gd", vim.lsp.buf.definition)
vim.keymap.set("n", "<leader>lr", require "nvchad.lsp.renamer")
vim.keymap.set("n", "gD", vim.lsp.buf.declaration)
-- vim.keymap.set("n", "<leader>wa", vim.lsp.buf.add_workspace_folder)
-- vim.keymap.set("n", "<leader>wr", vim.lsp.buf.remove_workspace_folder)
-- vim.keymap.set("n", "<leader>wl", function()
-- print(vim.inspect(vim.lsp.buf.list_workspace_folders()))
-- end)
-- vim.keymap.set("n", "<leader>D", vim.lsp.buf.type_definition)

--------------------- Mine / Prime ----------------------

-- File eXplorer -> doesn't work on nvchad?
-- vim.keymap.set("n", "<leader>x", vim.cmd.Ex)

-- Move commands within VISUAL mode
vim.keymap.set("v", "J", ":m '>+1<CR>gv=gv")
vim.keymap.set("v", "K", ":m '<-2<CR>gv=gv")

-- vim.api.nvim_set_keyvim.keymap.set("n", "<leader>tf", "<Plug>PlenaryTestFile", { noremap = false, silent = false })

vim.keymap.set("n", "J", "mzJ`z")
vim.keymap.set("n", "<C-d>", "<C-d>zz")
vim.keymap.set("n", "<C-u>", "<C-u>zz")
vim.keymap.set("n", "n", "nzzzv")
vim.keymap.set("n", "N", "Nzzzv")
vim.keymap.set("n", "=ap", "ma=ap'a")
vim.keymap.set("n", "<leader>zig", "<cmd>LspRestart<cr>")

-- greatest remap ever
vim.keymap.set("x", "<leader>p", [["_dP]])

-- Space y copies to SYSTEM clipboard
vim.keymap.set({ "n", "v" }, "<leader>y", [["+y]])
vim.keymap.set("n", "<leader>Y", [["+Y]])

-- Delete to void (i.e. preserving clipboard)
vim.keymap.set({ "n", "v" }, "<leader>d", '"_d')

-- Makes sure that vertical edits are performed
vim.keymap.set("i", "<C-c>", "<Esc>")

vim.keymap.set("n", "Q", "<nop>")
-- vim.keymap.set("n", "<C-f>", "<cmd>silent !tmux neww tmux-sessionizer<CR>")
vim.keymap.set("n", "<C-f>", function()
  require("conform").format { bufnr = 0 }
end)

vim.keymap.set("n", "<C-k>", "<cmd>cnext<CR>zz")
vim.keymap.set("n", "<C-j>", "<cmd>cprev<CR>zz")
-- vim.keymap.set("n", "<leader>k", "<cmd>lnext<CR>zz")
-- vim.keymap.set("n", "<leader>j", "<cmd>lprev<CR>zz")

vim.keymap.set("n", "<leader>r", [[:%s/\<<C-r><C-w>\>/<C-r><C-w>/gI<Left><Left><Left>]])
vim.keymap.set("n", "<leader>mx", "<cmd>!chmod +x %<CR>", { silent = true })

-- vim.keymap.set("n", "<leader>ee", "oif err != nil {<CR>}<Esc>Oreturn err<Esc>")
-- vim.keymap.set("n", "<leader>ea", 'oassert.NoError(err, "")<Esc>F";a')
-- vim.keymap.set("n", "<leader>ef", 'oif err != nil {<CR>}<Esc>Olog.Fatalf("error: %s\\n", err.Error())<Esc>jj')
-- vim.keymap.set("n", "<leader>el", 'oif err != nil {<CR>}<Esc>O.logger.Error("error", "error", err)<Esc>F.;i')

vim.keymap.set("n", "<leader><leader>", function()
  vim.cmd "so"
end)

vim.keymap.set("n", ";", ":", { desc = "CMD enter command mode" })
vim.keymap.set("i", "jk", "<ESC>")

-- telescope
vim.keymap.set("n", "<leader>s", "<cmd>Telescope live_grep<CR>", { desc = "telescope live grep" })
vim.keymap.set("n", "<leader>f", "<cmd>Telescope find_files<cr>", { desc = "telescope find files" })
vim.keymap.set("n", "<leader>fa", "<cmd>Telescope find_files follow=true no_ignore=true hidden=true<CR>")
vim.keymap.set("n", "<leader>g", "<cmd>Telescope git_files<CR>", { desc = "telescope git status" })

local builtin = require "telescope.builtin"
vim.keymap.set("n", "<leader>sw", function()
  local word = vim.fn.expand "<cword>"
  builtin.grep_string { search = word }
end)

-- vim.keymap.set("n", "<leader>fb", "<cmd>Telescope buffers<CR>", { desc = "telescope find buffers" })
-- vim.keymap.set("n", "<leader>fh", "<cmd>Telescope help_tags<CR>", { desc = "telescope help page" })
-- vim.keymap.set("n", "<leader>ma", "<cmd>Telescope marks<CR>", { desc = "telescope find marks" })
-- vim.keymap.set("n", "<leader>fo", "<cmd>Telescope oldfiles<CR>", { desc = "telescope find oldfiles" })
-- vim.keymap.set(
-- "n",
-- "<leader>fz",
-- "<cmd>Telescope current_buffer_fuzzy_find<CR>",
-- { desc = "telescope find in current buffer" }
-- )
-- vim.keymap.set("n", "<leader>cm", "<cmd>Telescope git_commits<CR>", { desc = "telescope git commits" })
-- vim.keymap.set("n", "<leader>gt", "<cmd>Telescope git_status<CR>", { desc = "telescope git status" })
-- vim.keymap.set("n", "<leader>pt", "<cmd>Telescope terms<CR>", { desc = "telescope pick hidden term" })
