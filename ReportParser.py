#!/usr/bin/python3

import json
import os
import bibtexparser
import bs4
import re

from bs4 import BeautifulSoup
from os import PathLike
from pathlib import Path, PosixPath
from string import capwords
from typing import Union
from unidecode import unidecode
from make_fmriprep_patterns import make_fmriprep_patterns
from parse_report import parse_report

# Compile regex patterns
# err_pat, getfig_pat, func_pat, fsub_pat, \
# title_pat, subtitle_pat, bp_pat = \
# tuple(map(re.compile, make_fmriprep_patterns().values()))
# 
# progdir = os.path.dirname(ReportParser.__file__)
# os.path.join(progdir, 'ReportParserPatterns.json')
# with  as jfile:
# # Artifacts to strip
# err_msg = ' '.join(('\s*Problem loading figure',
#                     '.*\.svg.',
#                     'If the link below works, please try reloading',
#                     'the report in your browser.\s*'))
# err_pat = re.compile(err_msg)
# getfig_pat = re.compile('\s*Get figure file\: .*.svg\s*')
# 
# # Headers
# func_pat = re.compile('Reports for: session .*, task .*\.')
# fsub_pat = re.compile('\n*Summary\n|\n*Confounds collected\n*')
# title_pat = re.compile('<div id=\"\w*\">\n')
# subtitle_pat = re.compile('.*\:\s{1}')
# 
# # Boilerplate section subtitles
# bp_pat = '|'.join(('\n*Anatomical data preprocessing\n',
#                    '\n*Functional data preprocessing\n',
#                    '\n*Copyright Waiver\n',
#                    '\n*References\n',
#                    '\n*Bibliography\n'))
# bp_pat = re.compile(bp_pat)
          




class ReportParser():
    def __init__(self,
                 src: Union[str, PathLike, PosixPath] = None,
                 json_file: Union[str, PathLike, PosixPath] = None,
                 # HTML parser engine
                 features: str = 'lxml',
                 # Carpet plot variables
                 cols: list = ['global_signal', 'csf',
                               'white_matter', 'csf_wm',
                               'framewise_displacement', 'std_dvars'],
                 **kwargs):
        """
        Instantiate a ``ReportParser`` object.
        """
        
        self.src = src
        self.features = features

        if kwargs is None:
            kwargs = {}
    # Compile regex patterns
        progdir = os.path.dirname(__file__)
        if json_file is None:
            json_file = os.path.join(progdir, 'ReportParserPatterns.json')
        with open(json_file, 'r') as jfile:
            patterns = json.load(jfile)
            jfile.close()
        __attributes__ = patterns    
        # __attributes__ = dict(err_pat=err_pat,
        #                       getfig_pat=getfig_pat,
        #                       title_pat=title_pat,
        #                       subtitle_pat=subtitle_pat,
        #                       func_pat=func_pat,
        #                       fsub_pat=fsub_pat,
        #                       bp_pat=bp_pat,
        #                       cols=cols)

        __attributes__.update(kwargs)
        self.__attributes__ = __attributes__
        
        [setattr(self, f'__{item[0]}__', item[1])
         for item in tuple(__attributes__.items())]

        self.__parse_report__ = parse_report


    def parse(self, src=None):
        
        self.src = src
        if self.src is None:
            parsed = {}
        else:
            parsed = self.__parse_report__(self.src, **self.__attributes__)

        return parsed

    @property
    def parsed(self):
        return self.__parse_report__(self.src)


    @property
    def sessions(self):
        return tuple(self.parsed['Functional'].keys())

    
    def __get_topdir__(self):
        psrc = Path(self.src)
        top = [p for p in psrc.parts
               if 'sub-' in p][0]
        ind = psrc.parts.index(top)
        return Path(os.path.join(*psrc.parts[:ind]))
    
    @property
    def topdir(self):
        """
        Returns the top BIDS directory.
        """
        
        return self.__get_topdir__()


    # Confounds and descriptive stats
    def __get_confounds__(self) -> None:
        filenames = [ses+'_desc-confounds_timeseries.tsv'
                     for ses in self.sessions]
        tables = [pd.read_csv(list(self.topdir.rglob('*'+n))[0],
                              sep='\t')
                  for n in filenames]
        return {self.sessions[table[0]]:table[1]
                for table in enumerate(tables)}

    @property
    def confounds(self):
        """
        Returns a dict containing confounds for each fMRI session.
        
        The dict is mapped using BIDS notation.
        Keys correspond to a specific session's confounds.
        Values are the respective confounds in a ``pandas.DataFrame``.
        For example:
            {'sub-<subject>_ses-<session>_task-<name>': confounds}
        """
        
        return self.__get_confounds__()


    def __get_confounds_summary__(self):
        return {self.sessions[conf[0]]:
                conf[1].loc[:, self.__cols__].describe()
                for conf in
                enumerate(tuple(self.confounds.values()))}

    @property
    def confounds_summary(self):
        """
        Returns descriptive stats displayed in the report's carpet plot.
        """

        return self.__get_confounds_summary__()


    def __get_func_summary__(self):
        fcp = self.parsed['Functional'].copy()
        [itm.pop('Confounds collected') for itm in tuple(fcp.values())]
        return pd.DataFrame((pd.Series(data=itm[1], name=itm[0])
                             for itm in tuple(tst['Functional'].items())))
    @property
    def func_summary(self):
        """
        Returns the report's "Functional" section as a dict.
        """

        return self.__get_func_summary__()


    def reset(self):
        """
        Returns the ``ReportParser`` instance with empty attributes.
        """

        d = vars(self)
        d['src'] = None
        d['parsed'] = {}
