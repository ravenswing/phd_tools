local options = {
  formatters_by_ft = {
    lua = { "stylua" },
    rust = { "rustfmt" },
    python = {
      -- To run the Ruff formatter.
      "ruff_format",
      -- To organize the imports.
      "ruff_organize_imports",
    },
    javascript = { "prettier" },
  },

  format_on_save = {
    -- These options will be passed to conform.format()
    timeout_ms = 500,
    lsp_fallback = true,
  },
}

return options
