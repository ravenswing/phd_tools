### Creates symmetry mates for each object loaded in PyMOL ###

Code works only for **PyMOL v1.7 and up**.
This chunk of code links F1 key to the function that generates symmetry mates within 5 Å from each object.
You can just add it to the `pymolrc.py` in your home folder.

```python
def create_symm_mates(selection):
# Creates symmetry mates within 5 A of the original structure.
# Colors the original in magenta, the rest in white.
# Linked to F1 key.
    for x in cmd.get_names(selection):
        cmd.show('lines', x)
        util.cbaw(x)
        cmd.show('cartoon', x)
        cmd.show('sticks', x + ' and hetatm')
        util.cbay(x + ' and hetatm')
        cmd.symexp(x + '_symm', x, x, 5)
        util.cbam(x)
cmd.set_key('F1', create_symm_mates, ['all'])
```