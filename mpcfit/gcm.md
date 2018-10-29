# Great Circle Motion

Raw code at this point.

Requires Python 3.6.  Venv strongly recommended if 3.6 or later is not your
distro standard Python.  Otherwise no dependencies outside of the Python
standard library.

It's a fairly literal port from some C code, so some parts not so Pythonic.

Needs work to become an easily usable package.

## Tests

Use the pytest package, but note there seem to be inconsistencies currently
between pytest documentation and what commands actually work.  One of
`py.test`, `pytest`, or `python -m pytest` might work.

Example, from the cloned gcm directory:

```
$ py.test --cov .
============================= test session starts ==============================
platform linux -- Python 3.6.0, pytest-3.2.3, py-1.4.34, pluggy-0.4.0
rootdir: /home/skeys/pp, inifile:
plugins: cov-2.5.1
collected 11 items                                                              

gcm_test.py ...........

----------- coverage: platform linux, python 3.6.0-final-0 -----------
Name          Stmts   Miss  Cover
---------------------------------
gcm.py          116      0   100%
gcm_test.py     100      0   100%
---------------------------------
TOTAL           216      0   100%


========================== 11 passed in 0.17 seconds ===========================
```

## Style

Passes flake8 with the .flake8 in the repo.
