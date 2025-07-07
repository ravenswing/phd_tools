vim.api.nvim_create_autocmd("BufEnter", {
  pattern = "*.py",
  callback = function()
    vim.fn.setreg("n", "GIdef main() -> None:<Enter>    ...<Enter>if __name__ == '__main__':<Enter>    main()")
  end,
})
