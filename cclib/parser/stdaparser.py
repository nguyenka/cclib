# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
#import collections

"""
    Parser for sTDA output.
    Include the following changes:
    ccio.py:
    from cclib.parser.stdaparser import sTDA
    (sTDA,      ["s  T  D  A"], True) to triggers
    parser/__init__.py:
    from cclib.parser.stdaparser import sTDA
"""

import numpy as np
from cclib.parser import logfileparser
from cclib.parser import utils

class sTDA(logfileparser.Logfile):
    """A sTDA log file."""

    def __init__(self, *args, **kwargs):
        super(sTDA, self).__init__(logname="sTDA", *args, **kwargs)

        # Flag for whether this calc is DFT.
        self.is_DFT = False
        self.singlet = False

    def __str__(self):
        """Return a string representation of the object."""
        return "sTDA output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'sTDA("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by Turbomole.

        The labels are standardized except for the first character being lowercase.
        """
        # TODO more work could be required, but we don't have any logfiles
        # with non-C1 symmetry.
        return label.capitalize()

    def before_parsing(self):

        self.periodic_table = utils.PeriodicTable()

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line[32:39] == "Version" in line:
            version = line.split()[2]
            self.metadata["package_version"] = version
        if "s T D A" in line :
            self.metadata['methods'].append("sTDA")
            self.metadata['success'] = False
        if "Welcome in nonlinear response" in line:
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            if "2PA" in line:
                self.metadata['methods'].append("QR-TPA")
        if 'atom   #          x             y             z' in line:
            atomcoords = []
            atomnos = []
            line = next(inputfile)
            while len(line.split()) == 6:
                #print(f'line {line}')
                atomnos.append(self.periodic_table.number[line.split()[0].capitalize()])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom")\
                        for x in line.split()[2:5]])
                line = next(inputfile)
            self.append_attribute('atomcoords', atomcoords)
            self.set_attribute('atomnos', atomnos)
            self.set_attribute('natom', len(atomcoords))
            #print(f'atomcoords: {atomcoords}')

        if "RKS-sTD-DFT" in line:
            self.singlet = True
            #print(f'self.singlet: {self.singlet}')

        if "roots found," in line:
            nstate = int(line.split()[0])
            #print(f'Number of excited states: {nstate}')
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            # Energy should be in cm-1...
            for n in range(nstate):
                line = next(inputfile)
                sym = 'A'
                if self.singlet:
                    mult='Singlet'
                else:
                    mult='???'
                symmetry = f"{mult}-{sym}"
                self.append_attribute("etsyms", symmetry)
                energy = utils.convertor(utils.float(line.split()[1]), 'eV', 'wavenumber')
                self.append_attribute("etenergies", energy)

                # Oscillator strength.
                oscillator_strength = utils.float(line.split()[3])
                #print(f'oscillator_strength: {oscillator_strength}')
                self.append_attribute("etoscs", oscillator_strength)
            
        if "#Delta ( " in line:
            #nstate = nstate + 1
            #print('found tpa info')
            energy = utils.convertor(utils.float(line.replace(')',' ').split()[2]), 'eV', 'wavenumber')
            print(f'energy {energy}')
            #self.append_attribute("etenergies", energy)
            #ab = next(inputfile) + next(inputfile) + next(inputfile)
            #print(f'ab {ab}')
            #ab_val = np.loadtxt(ab,usecols=(1,2,3), unpack=True)
            #print(f'ab_val {ab_val}')

        if "sTD-DFT done." in line:
            # End of module, set success flag.
            self.metadata['success'] = True
    def after_parsing(self):
        # atomcoords are parsed as a list of lists but it should be an array
        if hasattr(self, "atomcoords"):
            self.atomcoords = np.array(self.atomcoords)

