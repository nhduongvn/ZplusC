from contextlib import contextmanager

import ROOT


ROOT.gROOT.SetBatch(True)

@contextmanager
def open_root(filename='', mode=''):
    """
    A context manager for ROOT TFiles. The syntax mimics the built-in open function.
    """
    open_modes = {
        'r': 'READ',
        'w': 'RECREATE',
        'a': 'UPDATE',
        'r+': 'UPDATE',
        'w+': 'RECREATE',
        'a+': 'UPDATE',
    }

    if mode in open_modes:
        mode = open_modes[mode]
    else:
        mode = 'READ'

    root_file = ROOT.TFile.Open(filename, mode)
    if not root_file:
        raise IOError('Failed to open {}.'.format(filename))

    try:
        ROOT.SetOwnership(root_file, True)
        yield root_file
    finally:
        root_file.Close()

