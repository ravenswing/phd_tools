dofile(vim.g.base46_cache .. "telescope")

local actions = require "telescope.actions"

return {
  defaults = {
    selection_caret = " ",
    entry_prefix = " ",
    sorting_strategy = "ascending",
    layout_config = {
      horizontal = {
        prompt_position = "top",
        preview_width = 0.55,
      },
      width = 0.87,
      height = 0.80,
    },
    mappings = {
      n = { ["q"] = require("telescope.actions").close },
      i = {
        ["<C-h>"] = actions.select_horizontal,
        ["<esc>"] = actions.close,
      },
    },
  },

  extensions_list = { "themes", "terms" },
  extensions = {},
}
